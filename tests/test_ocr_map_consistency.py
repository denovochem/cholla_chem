"""
Test to ensure the shipped ocr_corrections_map.json matches what would be
generated from scratch. This catches cases where chemical_name_tokens.json
or substitution dicts are updated but the map wasn't regenerated.
"""

from __future__ import annotations

import json
from collections import deque
from pathlib import Path
from typing import Dict, List, Set

import pytest

from cholla_chem.utils.constants import (
    KEYBOARD_NEIGHBOR_SUBSTITUTIONS,
    OCR_SUBSTITUTIONS,
)

# Paths
DATAFILES_DIR = Path(__file__).resolve().parent.parent / "cholla_chem" / "datafiles"
CHEMICAL_NAME_TOKENS_PATH = DATAFILES_DIR / "chemical_name_tokens.json"
OCR_MAP_PATH = DATAFILES_DIR / "ocr_corrections_map.json"


def generate_substitution_dict(
    word: str, substitutions: Dict[str, List[str]], max_edits: int = 1
) -> Set[str]:
    """Generate all possible OCR/typo error variants of a word."""
    results = set()
    queue = deque([(word, 0)])
    visited = {word}

    while queue:
        current_s, depth = queue.popleft()
        if depth >= max_edits:
            continue

        for i in range(len(current_s)):
            for key, replacements in substitutions.items():
                key_len = len(key)
                if i + key_len > len(current_s):
                    continue

                if current_s[i : i + key_len] == key:
                    for replacement in replacements:
                        prefix = current_s[:i]
                        suffix = current_s[i + key_len :]
                        new_candidate = prefix + replacement + suffix

                        if new_candidate not in visited:
                            visited.add(new_candidate)
                            results.add(new_candidate)
                            queue.append((new_candidate, depth + 1))
    return results


def build_corrections_map_from_scratch(max_edits: int = 1) -> Dict[str, str]:
    """Regenerate the corrections map from source data."""
    with open(CHEMICAL_NAME_TOKENS_PATH, "r", encoding="utf-8") as f:
        chemical_name_tokens: List[str] = json.load(f)

    chemical_token_set = set(chemical_name_tokens)
    corrections_map: Dict[str, str] = {}

    substitution_sources = [
        OCR_SUBSTITUTIONS,
        KEYBOARD_NEIGHBOR_SUBSTITUTIONS,
    ]

    for source_dict in substitution_sources:
        for chemical_name_token in chemical_name_tokens:
            generated_errors = generate_substitution_dict(
                chemical_name_token, source_dict, max_edits=max_edits
            )

            for error in generated_errors:
                if len(error) <= 2:
                    continue
                if error in chemical_token_set:
                    continue

                clean_error_key = error.replace(" ", "")
                if clean_error_key not in corrections_map:
                    corrections_map[clean_error_key] = chemical_name_token

    return corrections_map


class TestOCRMapConsistency:
    """Ensure shipped OCR map is up to date."""

    def test_ocr_map_exists(self):
        """Verify the pre-computed map file exists."""
        assert OCR_MAP_PATH.exists(), (
            f"ocr_corrections_map.json not found at {OCR_MAP_PATH}. "
            "Run: python scripts/build_ocr_map.py"
        )

    def test_ocr_map_matches_generated(self):
        """Verify shipped map matches freshly generated map."""
        if not OCR_MAP_PATH.exists():
            pytest.skip("ocr_corrections_map.json does not exist")

        with open(OCR_MAP_PATH, "r", encoding="utf-8") as f:
            shipped_map = json.load(f)

        generated_map = build_corrections_map_from_scratch(max_edits=1)

        # Check for missing keys
        shipped_keys = set(shipped_map.keys())
        generated_keys = set(generated_map.keys())

        missing_in_shipped = generated_keys - shipped_keys
        extra_in_shipped = shipped_keys - generated_keys

        assert not missing_in_shipped, (
            f"ocr_corrections_map.json is missing {len(missing_in_shipped)} keys. "
            "Run: python scripts/build_ocr_map.py"
        )

        assert not extra_in_shipped, (
            f"ocr_corrections_map.json has {len(extra_in_shipped)} extra keys. "
            "Run: python scripts/build_ocr_map.py"
        )

        # Check for value mismatches
        mismatches = []
        for key in shipped_keys:
            if shipped_map[key] != generated_map[key]:
                mismatches.append(key)

        assert not mismatches, (
            f"ocr_corrections_map.json has {len(mismatches)} value mismatches. "
            "Run: python scripts/build_ocr_map.py"
        )
