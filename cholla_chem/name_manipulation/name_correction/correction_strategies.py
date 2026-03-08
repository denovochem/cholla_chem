from __future__ import annotations

import itertools
import json
import re
from abc import ABC, abstractmethod
from collections.abc import Iterator
from importlib import resources as importlib_resources
from typing import ClassVar, Dict, List, Optional, Tuple

from flashtext import KeywordProcessor

from cholla_chem.name_manipulation.name_correction.dataclasses import (
    Correction,
    CorrectionType,
    CorrectorConfig,
)
from cholla_chem.name_manipulation.name_correction.regexes import PATTERNS

# from cholla_chem.utils.logging_config import logger

FLASHTEXT_DICTS_PACKAGE = "cholla_chem.datafiles.flashtext_jsons"


class CorrectionStrategy(ABC):
    """
    Abstract base class for correction strategies.

    Subclass this to implement custom correction strategies.
    """

    @property
    @abstractmethod
    def name(self) -> str:
        """Return the name of this strategy."""
        pass

    @property
    @abstractmethod
    def correction_type(self) -> CorrectionType:
        """Return the type of corrections this strategy produces."""
        pass

    @abstractmethod
    def generate_candidates(
        self,
        text: str,
        num_corrections_applied: int = 0,
        config: Optional[CorrectorConfig] = None,
    ) -> Iterator[Tuple[str, List[Correction]]]:
        """
        Generate correction candidates from input text.

        Args:
            text: Input text to correct
            num_corrections_applied: Number of corrections already applied
            config: Corrector configuration

        Yields:
            Tuples of (corrected_text, list_of_corrections)
        """
        pass


class CharacterSubstitutionStrategy(CorrectionStrategy):
    """
    Strategy for correcting OCR character substitution errors using Aho-Corasick (FlashText).
    """

    @property
    def name(self) -> str:
        return "Character Substitution"

    @property
    def correction_type(self) -> CorrectionType:
        return CorrectionType.CHARACTER_SUBSTITUTION

    def __init__(
        self, max_edits: int = 1, substitutions: Optional[Dict[str, List[str]]] = None
    ):
        """
        Initialize with substitution map using FlashText for O(N) performance.

        Uses pre-computed corrections map for fast loading when no custom
        substitutions are provided.
        """
        self.keyword_processor = None
        self._substitutions = substitutions
        self._max_edits = max_edits

    def _initialize_keyword_processor(
        self,
        max_edits: int = 1,
    ):
        """
        Initialize the keyword processor with the given substitutions.
        """
        kp = KeywordProcessor()
        kp.non_word_boundaries = set()

        # Load pre-computed corrections map (fast path)
        with (
            importlib_resources.files(FLASHTEXT_DICTS_PACKAGE)
            .joinpath("substitutions_map.json")
            .open("r", encoding="utf-8") as f
        ):
            corrections_map = json.load(f)

        for error_key, correct_token in corrections_map.items():
            kp.add_keyword(error_key, correct_token)

        self.keyword_processor = kp

    def generate_candidates(
        self,
        text: str,
        num_corrections_applied: int = 0,
        config: Optional[CorrectorConfig] = None,
    ) -> Iterator[Tuple[str, List[Correction]]]:
        """
        Generate candidates by finding all matches and applying combinations of them.
        """
        if config is None:
            config = CorrectorConfig()

        if self.keyword_processor is None:
            self._initialize_keyword_processor(max_edits=self._max_edits)

        kp = self.keyword_processor
        assert kp is not None, "Keyword processor should be initialized"

        # 1. SQUASH & MAP
        clean_chars = []
        idx_map = []

        for i, char in enumerate(text):
            if not char.isspace():
                clean_chars.append(char)
                idx_map.append(i)

        clean_text = "".join(clean_chars)

        if not clean_text:
            return

        # 2. EXTRACT matches
        # matches = [('1,2-dichlorobenzene', 5, 12), ...]
        raw_matches = kp.extract_keywords(clean_text, span_info=True)

        if not raw_matches:
            return

        # 3. PREPARE Match Objects with Original Indices
        # We process these into a structured list so we can combine them easily
        mapped_matches = []

        for correct_token, start, end in raw_matches:
            # Map back to original text coordinates
            orig_start_index = idx_map[start]
            last_char_idx = idx_map[end - 1]
            orig_end_index = last_char_idx + 1

            original_match_string = text[orig_start_index:orig_end_index]

            match_data = {
                "replacement": correct_token,
                "start": orig_start_index,
                "end": orig_end_index,
                "original": original_match_string,
            }
            mapped_matches.append(match_data)

        # 4. GENERATE COMBINATIONS
        # By default, FlashText returns non-overlapping matches, so we don't need
        # complex overlap checks for the raw output.
        for r in range(
            1,
            min(
                len(mapped_matches),
                config.max_corrections_per_candidate - num_corrections_applied,
            )
            + 1,
        ):
            # itertools.combinations produces tuples of matches: (match1, match2)
            for combo in itertools.combinations(mapped_matches, r):
                # Sort combo by start index DESCENDING.
                # This is crucial: we must replace from end to start so
                # earlier indices don't shift when we modify the string.
                sorted_combo = sorted(combo, key=lambda x: x["start"], reverse=True)

                # Check for overlaps within this specific combination
                # (Paranoia check: FlashText usually handles this, but strictly speaking
                # if you have custom logic, this ensures safety)
                if self._check_overlap(sorted_combo):
                    continue

                # Apply corrections
                new_text_chars = list(text)  # Mutable char list
                corrections_list = []

                for match in sorted_combo:
                    # Replace in the char list
                    # We replace the slice with the new string
                    new_text_chars[match["start"] : match["end"]] = list(
                        match["replacement"]
                    )

                    # Create correction object
                    correction = Correction(
                        position=match["start"],
                        original=match["original"],
                        replacement=match["replacement"],
                        correction_type=self.correction_type,
                        description=f"OCR: '{match['original']}' → '{match['replacement']}'",
                    )
                    # Because we process in reverse for string building,
                    # we might want to insert at 0 to keep the list in reading order,
                    # or just append and reverse later.
                    corrections_list.append(correction)

                # Reassemble string
                corrected_text = "".join(new_text_chars)

                # Sort corrections back to reading order (optional, good for UI)
                corrections_list.sort(key=lambda x: x.position)

                yield corrected_text, corrections_list

    def _check_overlap(self, sorted_reverse_matches: List[Dict]) -> bool:
        """
        Check if any matches in the combination overlap.
        Input is sorted by start index DESCENDING.
        """
        # range (start, end)
        # B is before A in the list (because reverse sort), so B has higher index
        for i in range(len(sorted_reverse_matches) - 1):
            later_match = sorted_reverse_matches[i]
            earlier_match = sorted_reverse_matches[i + 1]

            # If the earlier match ends after the later match starts -> Overlap
            if earlier_match["end"] > later_match["start"]:
                return True
        return False

    @staticmethod
    def _positions_overlap(points: List[Tuple[int, str, List[str]]]) -> bool:
        """Check if any substitution positions overlap."""
        ranges = [(p[0], p[0] + len(p[1])) for p in points]
        ranges.sort()
        for i in range(len(ranges) - 1):
            if ranges[i][1] > ranges[i + 1][0]:
                return True
        return False


class CharacterInsertionStrategy(CorrectionStrategy):
    @property
    def name(self) -> str:
        return "Character Insertion"

    @property
    def correction_type(self) -> CorrectionType:
        return CorrectionType.CHARACTER_INSERTION

    def __init__(self, max_edits: int = 1):
        self.keyword_processor = None
        self._max_edits = max_edits

    def _initialize_keyword_processor(self, max_edits: int = 1):
        kp = KeywordProcessor()
        kp.non_word_boundaries = set()

        # Load pre-computed corrections map (fast path)
        with (
            importlib_resources.files(FLASHTEXT_DICTS_PACKAGE)
            .joinpath("insertions_map.json")
            .open("r", encoding="utf-8") as f
        ):
            corrections_map = json.load(f)

        for error_key, correct_token in corrections_map.items():
            kp.add_keyword(error_key, correct_token)

        self.keyword_processor = kp

    def generate_candidates(
        self,
        text: str,
        num_corrections_applied: int = 0,
        config: Optional[CorrectorConfig] = None,
    ) -> Iterator[Tuple[str, List[Correction]]]:
        if config is None:
            config = CorrectorConfig()

        if self.keyword_processor is None:
            self._initialize_keyword_processor(max_edits=self._max_edits)

        kp = self.keyword_processor
        assert kp is not None, "Keyword processor should be initialized"

        clean_chars = []
        idx_map = []

        for i, char in enumerate(text):
            if not char.isspace():
                clean_chars.append(char)
                idx_map.append(i)

        clean_text = "".join(clean_chars)
        if not clean_text:
            return

        raw_matches = kp.extract_keywords(clean_text, span_info=True)
        if not raw_matches:
            return

        mapped_matches = []
        for correct_token, start, end in raw_matches:
            orig_start_index = idx_map[start]
            last_char_idx = idx_map[end - 1]
            orig_end_index = last_char_idx + 1

            original_match_string = text[orig_start_index:orig_end_index]
            if original_match_string == correct_token:
                continue

            mapped_matches.append(
                {
                    "replacement": correct_token,
                    "start": orig_start_index,
                    "end": orig_end_index,
                    "original": original_match_string,
                }
            )

        for r in range(
            1,
            min(
                len(mapped_matches),
                config.max_corrections_per_candidate - num_corrections_applied,
            )
            + 1,
        ):
            for combo in itertools.combinations(mapped_matches, r):
                sorted_combo = sorted(combo, key=lambda x: x["start"], reverse=True)
                if self._check_overlap(sorted_combo):
                    continue

                new_text_chars = list(text)
                corrections_list: List[Correction] = []

                for match in sorted_combo:
                    new_text_chars[match["start"] : match["end"]] = list(
                        match["replacement"]
                    )
                    corrections_list.append(
                        Correction(
                            position=match["start"],
                            original=match["original"],
                            replacement=match["replacement"],
                            correction_type=self.correction_type,
                            description=f"Insertion: '{match['original']}' → '{match['replacement']}'",
                        )
                    )

                corrected_text = "".join(new_text_chars)
                corrections_list.sort(key=lambda x: x.position)
                yield corrected_text, corrections_list

    def _check_overlap(self, sorted_reverse_matches: List[Dict]) -> bool:
        for i in range(len(sorted_reverse_matches) - 1):
            later_match = sorted_reverse_matches[i]
            earlier_match = sorted_reverse_matches[i + 1]
            if earlier_match["end"] > later_match["start"]:
                return True
        return False


class CharacterDeletionStrategy(CorrectionStrategy):
    @property
    def name(self) -> str:
        return "Character Deletion"

    @property
    def correction_type(self) -> CorrectionType:
        return CorrectionType.CHARACTER_DELETION

    def __init__(self, max_edits: int = 1):
        self.keyword_processor = None
        self._max_edits = max_edits

    def _initialize_keyword_processor(self, max_edits: int = 1):
        kp = KeywordProcessor()
        kp.non_word_boundaries = set()

        # Load pre-computed corrections map (fast path)
        with (
            importlib_resources.files(FLASHTEXT_DICTS_PACKAGE)
            .joinpath("deletions_map.json")
            .open("r", encoding="utf-8") as f
        ):
            corrections_map = json.load(f)

        for error_key, correct_token in corrections_map.items():
            kp.add_keyword(error_key, correct_token)

        self.keyword_processor = kp

    def generate_candidates(
        self,
        text: str,
        num_corrections_applied: int = 0,
        config: Optional[CorrectorConfig] = None,
    ) -> Iterator[Tuple[str, List[Correction]]]:
        if config is None:
            config = CorrectorConfig()

        if self.keyword_processor is None:
            self._initialize_keyword_processor(max_edits=self._max_edits)

        kp = self.keyword_processor
        assert kp is not None, "Keyword processor should be initialized"

        clean_chars = []
        idx_map = []

        for i, char in enumerate(text):
            if not char.isspace():
                clean_chars.append(char)
                idx_map.append(i)

        clean_text = "".join(clean_chars)
        if not clean_text:
            return

        raw_matches = kp.extract_keywords(clean_text, span_info=True)
        if not raw_matches:
            return

        mapped_matches = []
        for correct_token, start, end in raw_matches:
            orig_start_index = idx_map[start]
            last_char_idx = idx_map[end - 1]
            orig_end_index = last_char_idx + 1

            original_match_string = text[orig_start_index:orig_end_index]
            if original_match_string == correct_token:
                continue

            mapped_matches.append(
                {
                    "replacement": correct_token,
                    "start": orig_start_index,
                    "end": orig_end_index,
                    "original": original_match_string,
                }
            )

        for r in range(
            1,
            min(
                len(mapped_matches),
                config.max_corrections_per_candidate - num_corrections_applied,
            )
            + 1,
        ):
            for combo in itertools.combinations(mapped_matches, r):
                sorted_combo = sorted(combo, key=lambda x: x["start"], reverse=True)
                if self._check_overlap(sorted_combo):
                    continue

                new_text_chars = list(text)
                corrections_list: List[Correction] = []

                for match in sorted_combo:
                    new_text_chars[match["start"] : match["end"]] = list(
                        match["replacement"]
                    )
                    corrections_list.append(
                        Correction(
                            position=match["start"],
                            original=match["original"],
                            replacement=match["replacement"],
                            correction_type=self.correction_type,
                            description=f"Deletion: '{match['original']}' → '{match['replacement']}'",
                        )
                    )

                corrected_text = "".join(new_text_chars)
                corrections_list.sort(key=lambda x: x.position)
                yield corrected_text, corrections_list

    def _check_overlap(self, sorted_reverse_matches: List[Dict]) -> bool:
        for i in range(len(sorted_reverse_matches) - 1):
            later_match = sorted_reverse_matches[i]
            earlier_match = sorted_reverse_matches[i + 1]
            if earlier_match["end"] > later_match["start"]:
                return True
        return False


class CharacterTranspositionStrategy(CorrectionStrategy):
    @property
    def name(self) -> str:
        return "Character Transposition"

    @property
    def correction_type(self) -> CorrectionType:
        return CorrectionType.CHARACTER_TRANSPOSITION

    def __init__(self, max_edits: int = 1):
        self.keyword_processor = None
        self._max_edits = max_edits

    def _initialize_keyword_processor(self, max_edits: int = 1):
        kp = KeywordProcessor()
        kp.non_word_boundaries = set()

        # Load pre-computed corrections map (fast path)
        with (
            importlib_resources.files(FLASHTEXT_DICTS_PACKAGE)
            .joinpath("transpositions_map.json")
            .open("r", encoding="utf-8") as f
        ):
            corrections_map = json.load(f)

        for error_key, correct_token in corrections_map.items():
            kp.add_keyword(error_key, correct_token)

        self.keyword_processor = kp

    def generate_candidates(
        self,
        text: str,
        num_corrections_applied: int = 0,
        config: Optional[CorrectorConfig] = None,
    ) -> Iterator[Tuple[str, List[Correction]]]:
        if config is None:
            config = CorrectorConfig()

        if self.keyword_processor is None:
            self._initialize_keyword_processor(max_edits=self._max_edits)

        kp = self.keyword_processor
        assert kp is not None, "Keyword processor should be initialized"

        clean_chars = []
        idx_map = []

        for i, char in enumerate(text):
            if not char.isspace():
                clean_chars.append(char)
                idx_map.append(i)

        clean_text = "".join(clean_chars)
        if not clean_text:
            return

        raw_matches = kp.extract_keywords(clean_text, span_info=True)
        if not raw_matches:
            return

        mapped_matches = []
        for correct_token, start, end in raw_matches:
            orig_start_index = idx_map[start]
            last_char_idx = idx_map[end - 1]
            orig_end_index = last_char_idx + 1

            original_match_string = text[orig_start_index:orig_end_index]
            if original_match_string == correct_token:
                continue

            mapped_matches.append(
                {
                    "replacement": correct_token,
                    "start": orig_start_index,
                    "end": orig_end_index,
                    "original": original_match_string,
                }
            )

        for r in range(
            1, min(len(mapped_matches), config.max_corrections_per_candidate) + 1
        ):
            for combo in itertools.combinations(mapped_matches, r):
                sorted_combo = sorted(combo, key=lambda x: x["start"], reverse=True)
                if self._check_overlap(sorted_combo):
                    continue

                new_text_chars = list(text)
                corrections_list: List[Correction] = []

                for match in sorted_combo:
                    new_text_chars[match["start"] : match["end"]] = list(
                        match["replacement"]
                    )
                    corrections_list.append(
                        Correction(
                            position=match["start"],
                            original=match["original"],
                            replacement=match["replacement"],
                            correction_type=self.correction_type,
                            description=f"Transposition: '{match['original']}' → '{match['replacement']}'",
                        )
                    )

                corrected_text = "".join(new_text_chars)
                corrections_list.sort(key=lambda x: x.position)
                yield corrected_text, corrections_list

    def _check_overlap(self, sorted_reverse_matches: List[Dict]) -> bool:
        for i in range(len(sorted_reverse_matches) - 1):
            later_match = sorted_reverse_matches[i]
            earlier_match = sorted_reverse_matches[i + 1]
            if earlier_match["end"] > later_match["start"]:
                return True
        return False


class PunctuationRestorationStrategy(CorrectionStrategy):
    """
    Strategy for restoring missing punctuation in chemical names.

    Handles:
    - Missing hyphens between locants and substituent names
    - Missing commas between locants
    - Missing brackets around stereochemistry
    """

    @property
    def name(self) -> str:
        return "Punctuation Restoration"

    @property
    def correction_type(self) -> CorrectionType:
        return CorrectionType.PUNCTUATION_RESTORATION

    # Patterns for missing punctuation
    MISSING_HYPHEN_PATTERNS: ClassVar[List[Tuple[re.Pattern, str, str]]] = [
        # digit followed directly by letter (missing hyphen after locant)
        (re.compile(r"(\d)([a-zA-Z])"), r"\1-\2", "digit-letter"),
        # letter followed by digit (missing hyphen before locant)
        (re.compile(r"([a-zA-Z])(\d)"), r"\1-\2", "letter-digit"),
    ]

    MISSING_COMMA_PATTERNS: ClassVar[List[Tuple[re.Pattern, str, str]]] = [
        # Adjacent digits that should be comma-separated locants
        (re.compile(r"(\d)\s+(\d)"), r"\1,\2", "spaced-digits"),
    ]

    def generate_candidates(
        self,
        text: str,
        num_corrections_applied: int = 0,
        config: Optional[CorrectorConfig] = None,
    ) -> Iterator[Tuple[str, List[Correction]]]:
        """Generate candidates by restoring missing punctuation."""
        yield text, []  # Original

        # Try hyphen restoration
        for pattern, replacement, desc in self.MISSING_HYPHEN_PATTERNS:
            matches = list(pattern.finditer(text))

            for match in matches:
                new_text = (
                    text[: match.start()]
                    + pattern.sub(replacement, match.group())
                    + text[match.end() :]
                )

                if new_text != text:
                    correction = Correction(
                        position=match.start(),
                        original=match.group(),
                        replacement=pattern.sub(replacement, match.group()),
                        correction_type=self.correction_type,
                        description=f"Added hyphen ({desc})",
                    )
                    yield new_text, [correction]

        # Try comma restoration
        for pattern, replacement, desc in self.MISSING_COMMA_PATTERNS:
            matches = list(pattern.finditer(text))

            for match in matches:
                new_text = (
                    text[: match.start()]
                    + pattern.sub(replacement, match.group())
                    + text[match.end() :]
                )

                if new_text != text:
                    correction = Correction(
                        position=match.start(),
                        original=match.group(),
                        replacement=pattern.sub(replacement, match.group()),
                        correction_type=self.correction_type,
                        description=f"Added comma ({desc})",
                    )
                    yield new_text, [correction]


class BracketBalancingStrategy(CorrectionStrategy):
    """
    Strategy for balancing brackets in chemical names.

    Attempts to fix:
    - Missing closing brackets
    - Missing opening brackets
    - Mismatched bracket types
    """

    @property
    def name(self) -> str:
        return "Bracket Balancing"

    @property
    def correction_type(self) -> CorrectionType:
        return CorrectionType.BRACKET_BALANCING

    BRACKET_PAIRS: ClassVar[List[Tuple[str, str]]] = [
        ("(", ")"),
        ("[", "]"),
        ("{", "}"),
    ]

    def generate_candidates(
        self,
        text: str,
        num_corrections_applied: int = 0,
        config: Optional[CorrectorConfig] = None,
    ) -> Iterator[Tuple[str, List[Correction]]]:
        """Generate candidates by balancing brackets."""
        yield text, []  # Original

        for open_b, close_b in self.BRACKET_PAIRS:
            open_count = text.count(open_b)
            close_count = text.count(close_b)

            if open_count > close_count:
                # Missing closing brackets - add at end
                diff = open_count - close_count
                new_text = text + close_b * diff
                correction = Correction(
                    position=len(text),
                    original="",
                    replacement=close_b * diff,
                    correction_type=self.correction_type,
                    description=f"Added {diff} closing '{close_b}'",
                )
                yield new_text, [correction]

            elif close_count > open_count:
                # Missing opening brackets - add at start
                diff = close_count - open_count
                new_text = open_b * diff + text
                correction = Correction(
                    position=0,
                    original="",
                    replacement=open_b * diff,
                    correction_type=self.correction_type,
                    description=f"Added {diff} opening '{open_b}'",
                )
                yield new_text, [correction]


class LocantCorrectionStrategy(CorrectionStrategy):
    """
    Regex-based strategy for correcting chemical locants, stereochemistry,
    and punctuation errors (l->1, .->, missing commas).
    """

    @property
    def name(self) -> str:
        return "Locant Regex Correction"

    @property
    def correction_type(self) -> CorrectionType:
        return CorrectionType.LOCANT_FIX

    def generate_candidates(
        self,
        text: str,
        num_corrections_applied: int = 0,
        config: Optional[CorrectorConfig] = None,
    ) -> Iterator[Tuple[str, List[Correction]]]:
        # We apply these patterns sequentially to the string.
        # Order matters: Fix characters -> Fix punctuation -> Fix formatting.

        # In this implementation, we apply the corrections sequentially
        # to build one final "Clean" string, yielding each step.

        current_text = text

        for pattern, replacement, desc in PATTERNS:
            # We use a loop to catch all instances of one error type (e.g. multiple 'l's)
            # before moving to the next pattern type.

            matches = list(re.finditer(pattern, current_text))
            if not matches:
                continue

            corrections_batch = []

            # We must apply replacements from right to left (reverse)
            # so indices don't shift for the remaining matches in this batch.
            for match in reversed(matches):
                start, end = match.span()
                original_str = match.group()

                # specific handling for regex group replacement (like \1,\2)
                if isinstance(replacement, str) and "\\" in replacement:
                    new_str = match.expand(replacement)
                else:
                    new_str = replacement

                # Create Correction Object
                correction = Correction(
                    position=start,
                    original=original_str,
                    replacement=new_str,
                    correction_type=self.correction_type,
                    description=desc,
                )
                corrections_batch.append(correction)

                # Mutate the string
                current_text = current_text[:start] + new_str + current_text[end:]

            # Since we processed reverse, reverse back for the UI/List
            corrections_batch.reverse()

            yield current_text, corrections_batch
