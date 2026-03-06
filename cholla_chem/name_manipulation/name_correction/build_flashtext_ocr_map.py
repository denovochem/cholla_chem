from __future__ import annotations

import json
import os
from collections import deque
from pathlib import Path
from typing import Dict, List, Set

# Import substitution dicts from the package
from cholla_chem.utils.constants import (
    KEYBOARD_NEIGHBOR_SUBSTITUTIONS,
    OCR_SUBSTITUTIONS,
)

BASE_DIR = Path(__file__).resolve().parent
CHEMICAL_NAME_TOKENS_PATH = os.path.abspath(
    BASE_DIR.parent.parent / "datafiles" / "chemical_name_tokens.json"
)


def generate_substitution_dict(
    word: str, substitutions: Dict[str, List[str]], max_edits: int = 1
) -> Set[str]:
    """
    Generate all possible substitution error variants of a word.

    Args:
        word: The word to generate substitution errors for
        substitutions: Dictionary of character substitutions
        max_edits: Maximum number of substitutions to generate

    Returns:
        Set of all possible substitution error variants
    """
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


def generate_insertion_errors(word: str, max_edits: int = 1) -> Set[str]:
    """
    Generate all possible insertion error variants of a word.

    Args:
        word: The word to generate insertion errors for
        max_edits: Maximum number of insertions to generate

    Returns:
        Set of all possible insertion error variants
    """
    results: Set[str] = set()
    queue = deque([(word, 0)])
    visited = {word}

    while queue:
        current_s, depth = queue.popleft()
        if depth >= max_edits:
            continue

        for i in range(len(current_s) + 1):
            left = current_s[i - 1] if i > 0 else None
            right = current_s[i] if i < len(current_s) else None

            insertion_candidates: Set[str] = set()

            if left is not None:
                insertion_candidates.add(left)
                insertion_candidates.update(
                    KEYBOARD_NEIGHBOR_SUBSTITUTIONS.get(left, [])
                )

            if right is not None:
                insertion_candidates.add(right)
                insertion_candidates.update(
                    KEYBOARD_NEIGHBOR_SUBSTITUTIONS.get(right, [])
                )

            for ch in insertion_candidates:
                new_candidate = current_s[:i] + ch + current_s[i:]
                if new_candidate not in visited:
                    visited.add(new_candidate)
                    results.add(new_candidate)
                    queue.append((new_candidate, depth + 1))

    return results


def generate_deletion_errors(word: str, max_edits: int = 1) -> Set[str]:
    """
    Generate all possible deletion error variants of a word.

    Args:
        word: The word to generate deletion errors for
        max_edits: Maximum number of deletions to generate

    Returns:
        Set of all possible deletion error variants
    """
    results: Set[str] = set()
    queue = deque([(word, 0)])
    visited = {word}

    while queue:
        current_s, depth = queue.popleft()
        if depth >= max_edits:
            continue

        for i in range(len(current_s)):
            new_candidate = current_s[:i] + current_s[i + 1 :]
            if new_candidate not in visited:
                visited.add(new_candidate)
                results.add(new_candidate)
                queue.append((new_candidate, depth + 1))

    return results


def generate_transposition_errors(word: str, max_edits: int = 1) -> Set[str]:
    """
    Generate all possible transposition error variants of a word.

    Args:
        word: The word to generate transposition errors for
        max_edits: Maximum number of transpositions to generate

    Returns:
        Set of all possible transposition error variants
    """
    results: Set[str] = set()
    queue = deque([(word, 0)])
    visited = {word}

    while queue:
        current_s, depth = queue.popleft()
        if depth >= max_edits:
            continue

        # swap adjacent characters
        for i in range(len(current_s) - 1):
            if current_s[i] == current_s[i + 1]:
                continue
            swapped = (
                current_s[:i] + current_s[i + 1] + current_s[i] + current_s[i + 2 :]
            )

            if swapped not in visited:
                visited.add(swapped)
                results.add(swapped)
                queue.append((swapped, depth + 1))

    return results


def build_substitution_corrections_map(max_edits: int = 1) -> Dict[str, str]:
    """
    Build the complete error -> correct_token mapping.

    Returns:
        Dict mapping error strings (squashed) to correct chemical tokens.
    """
    with open(CHEMICAL_NAME_TOKENS_PATH, "r", encoding="utf-8") as f:
        chemical_name_tokens: List[str] = json.load(f)

    chemical_token_set = set(chemical_name_tokens)
    corrections_map: Dict[str, str] = {}

    substitution_sources = [
        ("OCR", OCR_SUBSTITUTIONS),
        ("Typo", KEYBOARD_NEIGHBOR_SUBSTITUTIONS),
    ]

    for source_name, source_dict in substitution_sources:
        count = 0
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
                # Only add if not already present (first mapping wins)
                if clean_error_key not in corrections_map:
                    corrections_map[clean_error_key] = chemical_name_token
                    count += 1

        print(f"[{source_name}] Added {count} new mappings")

    print(f"Total mappings: {len(corrections_map)}")
    return corrections_map


def build_insertion_corrections_map(max_edits: int = 1) -> Dict[str, str]:
    """
    Build insertion error -> correct_token mapping.
    """
    with open(CHEMICAL_NAME_TOKENS_PATH, "r", encoding="utf-8") as f:
        chemical_name_tokens: List[str] = json.load(f)

    chemical_token_set = set(chemical_name_tokens)
    corrections_map: Dict[str, str] = {}

    count = 0
    for token in chemical_name_tokens:
        generated_errors = generate_insertion_errors(token, max_edits=max_edits)

        for error in generated_errors:
            if len(error) <= 2:
                continue
            if error in chemical_token_set:
                continue

            clean_error_key = error.replace(" ", "")
            if clean_error_key not in corrections_map:
                corrections_map[clean_error_key] = token
                count += 1

    print(f"[Insertion] Added {count} new mappings")
    print(f"[Insertion] Total mappings: {len(corrections_map)}")
    return corrections_map


def build_deletion_corrections_map(max_edits: int = 1) -> Dict[str, str]:
    """
    Build deletion error -> correct_token mapping.
    """
    with open(CHEMICAL_NAME_TOKENS_PATH, "r", encoding="utf-8") as f:
        chemical_name_tokens: List[str] = json.load(f)

    chemical_token_set = set(chemical_name_tokens)
    corrections_map: Dict[str, str] = {}

    count = 0
    for token in chemical_name_tokens:
        generated_errors = generate_deletion_errors(token, max_edits=max_edits)

        for error in generated_errors:
            if len(error) <= 2:
                continue
            if error in chemical_token_set:
                continue

            clean_error_key = error.replace(" ", "")
            if clean_error_key not in corrections_map:
                corrections_map[clean_error_key] = token
                count += 1

    print(f"[Deletion] Added {count} new mappings")
    print(f"[Deletion] Total mappings: {len(corrections_map)}")
    return corrections_map


def build_transposition_corrections_map(max_edits: int = 1) -> Dict[str, str]:
    """
    Build transposition error -> correct_token mapping.
    """
    with open(CHEMICAL_NAME_TOKENS_PATH, "r", encoding="utf-8") as f:
        chemical_name_tokens: List[str] = json.load(f)

    chemical_token_set = set(chemical_name_tokens)
    corrections_map: Dict[str, str] = {}

    count = 0
    for token in chemical_name_tokens:
        generated_errors = generate_transposition_errors(token, max_edits=max_edits)

        for error in generated_errors:
            if len(error) <= 2:
                continue
            if error in chemical_token_set:
                continue

            clean_error_key = error.replace(" ", "")
            if clean_error_key not in corrections_map:
                corrections_map[clean_error_key] = token
                count += 1

    print(f"[Transposition] Added {count} new mappings")
    print(f"[Transposition] Total mappings: {len(corrections_map)}")
    return corrections_map
