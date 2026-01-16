from __future__ import annotations
import re
from typing import ClassVar, Dict, Optional

from placeholder_name.name_manipulation.name_correction.dataclasses import (
    CorrectionCandidate,
    CorrectorConfig,
)


class ChemicalNameScorer:
    """
    Scores chemical name candidates based on various heuristics.

    The scoring system uses multiple factors:
    - Bracket balance and validity
    - Presence of known chemical prefixes/suffixes
    - Valid locant patterns
    - Stereochemistry notation validity
    - Number of corrections (fewer is better)
    - Character distribution plausibility
    """

    # Weights for different scoring components
    WEIGHTS: ClassVar[Dict[str, float]] = {
        "bracket_balance": 0.20,
        "locant_validity": 0.15,
        "stereochemistry_validity": 0.10,
        "correction_penalty": 0.20,
        "character_plausibility": 0.10,
        "length_similarity": 0.25,
    }

    # Precompiled regex patterns
    LOCANT_PATTERN: ClassVar[re.Pattern] = re.compile(r"^\d+[,\d]*-|,\d+-|-\d+,")

    STEREOCHEM_PATTERN: ClassVar[re.Pattern] = re.compile(
        r"^\(R\)-|^\(S\)-|^\(E\)-|^\(Z\)-|^D-|^L-|^DL-|^rac-|^rel-|"
        r"\(R\)|\(S\)|\(E\)|\(Z\)|"
        r"-\(R\)-|-\(S\)-|-\(E\)-|-\(Z\)-"
    )

    VALID_CHARS_PATTERN: ClassVar[re.Pattern] = re.compile(
        r"^[a-zA-Z0-9\-\(\)\[\]\{\},\.\'\"\s]+$"
    )

    def __init__(self, config: Optional[CorrectorConfig] = None):
        """
        Initialize the scorer.

        Args:
            config: Optional corrector configuration
        """
        self.config = config or CorrectorConfig()

    def score(self, candidate: CorrectionCandidate) -> CorrectionCandidate:
        """
        Calculate the composite score for a candidate.

        Args:
            candidate: The candidate to score

        Returns:
            The same candidate with updated score and score_components
        """
        components: Dict[str, float] = {}

        # Calculate individual components
        components["bracket_balance"] = self._score_bracket_balance(candidate.name)
        components["locant_validity"] = self._score_locants(candidate.name)
        components["stereochemistry_validity"] = self._score_stereochemistry(
            candidate.name
        )
        components["correction_penalty"] = self._score_correction_count(
            candidate.num_corrections
        )
        components["character_plausibility"] = self._score_characters(candidate.name)
        components["length_similarity"] = self._score_length_similarity(
            candidate.name, candidate.original_name
        )

        # Calculate weighted sum
        total_score = sum(components[key] * self.WEIGHTS[key] for key in self.WEIGHTS)

        candidate.score = total_score
        candidate.score_components = components

        return candidate

    def _score_bracket_balance(self, name: str) -> float:
        """
        Score based on bracket balance.

        Returns 1.0 if all brackets are balanced, 0.0 if severely unbalanced.
        """
        bracket_pairs = [("(", ")"), ("[", "]"), ("{", "}")]
        total_imbalance = 0

        for open_b, close_b in bracket_pairs:
            open_count = name.count(open_b)
            close_count = name.count(close_b)
            total_imbalance += abs(open_count - close_count)

        # Also check for proper nesting
        if not self._check_bracket_nesting(name):
            total_imbalance += 2

        # Convert imbalance to score (0-1)
        return max(0.0, 1.0 - (total_imbalance * 0.25))

    def _check_bracket_nesting(self, name: str) -> bool:
        """Check if brackets are properly nested."""
        stack = []
        bracket_map = {")": "(", "]": "[", "}": "{"}

        for char in name:
            if char in "([{":
                stack.append(char)
            elif char in ")]}":
                if not stack or stack[-1] != bracket_map[char]:
                    return False
                stack.pop()

        return len(stack) == 0

    def _score_locants(self, name: str) -> float:
        """Score based on valid locant patterns."""
        # Check for proper locant formatting (e.g., "2,3-" or "1-")
        locant_matches = re.findall(r"\d+[,\d]*-", name)

        if not locant_matches:
            return 0.5  # Neutral score if no locants

        score = 0.7

        for match in locant_matches:
            # Check if locant numbers are reasonable (1-99 typically)
            numbers = re.findall(r"\d+", match)
            for num_str in numbers:
                num = int(num_str)
                if 1 <= num <= 50:
                    score += 0.05
                elif num > 100:
                    score -= 0.1

        return max(0.0, min(1.0, score))

    def _score_stereochemistry(self, name: str) -> float:
        """Score based on stereochemistry notation validity."""
        # Look for stereochemistry patterns
        stereo_matches = self.STEREOCHEM_PATTERN.findall(name)

        if not stereo_matches:
            return 0.5  # Neutral if no stereochemistry

        # Valid stereochemistry patterns boost score
        return min(1.0, 0.7 + len(stereo_matches) * 0.1)

    def _score_correction_count(self, num_corrections: int) -> float:
        """
        Score based on number of corrections (fewer is better).

        Follows the principle of minimal correction.
        """
        if num_corrections == 0:
            return 1.0
        elif num_corrections == 1:
            return 0.9
        elif num_corrections == 2:
            return 0.7
        elif num_corrections == 3:
            return 0.5
        else:
            return max(0.1, 1.0 - num_corrections * 0.15)

    def _score_characters(self, name: str) -> float:
        """Score based on character plausibility for chemical names."""
        if not name:
            return 0.0

        # Check for valid character set
        if not self.VALID_CHARS_PATTERN.match(name):
            return 0.3

        # Calculate letter/digit ratio (chemical names are mostly letters)
        letters = sum(1 for c in name if c.isalpha())
        total = len(name)

        if total == 0:
            return 0.0

        letter_ratio = letters / total

        # Ideal ratio is mostly letters with some digits
        if 0.6 <= letter_ratio <= 0.95:
            return 0.9
        elif 0.4 <= letter_ratio < 0.6:
            return 0.6
        elif letter_ratio > 0.95:
            return 0.8
        else:
            return 0.4

    def _score_length_similarity(self, corrected: str, original: str) -> float:
        """Score based on length similarity to original."""
        if not original:
            return 0.5

        len_diff = abs(len(corrected) - len(original))
        max_len = max(len(corrected), len(original))

        if max_len == 0:
            return 1.0

        similarity = 1.0 - (len_diff / max_len)
        return max(0.0, similarity)
