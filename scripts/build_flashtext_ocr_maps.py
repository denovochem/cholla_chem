#!/usr/bin/env python
"""
Build script to pre-compute OCR correction mappings.

Run this whenever chemical_name_tokens.json or substitution dicts change:
    python scripts/build_ocr_map.py
"""

import json
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent

sys.path.insert(0, str(PROJECT_ROOT))

from cholla_chem.name_manipulation.name_correction.build_flashtext_ocr_map import (  # noqa: E402
    build_deletion_corrections_map,
    build_insertion_corrections_map,
    build_substitution_corrections_map,
    build_transposition_corrections_map,
)

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
DATAFILES_DIR = PROJECT_ROOT / "cholla_chem" / "datafiles" / "flashtext_jsons"
SUBSTITUTION_OUTPUT_PATH = DATAFILES_DIR / "substitutions_map.json"
DELETION_OUTPUT_PATH = DATAFILES_DIR / "deletions_map.json"
INSERTION_OUTPUT_PATH = DATAFILES_DIR / "insertions_map.json"
TRANSPOSITION_OUTPUT_PATH = DATAFILES_DIR / "transpositions_map.json"


def main():
    substitution_corrections_map = build_substitution_corrections_map(max_edits=1)
    with open(SUBSTITUTION_OUTPUT_PATH, "w", encoding="utf-8") as f:
        json.dump(substitution_corrections_map, f, ensure_ascii=False, indent=2)

    print(
        f"\nWrote {len(substitution_corrections_map)} mappings to {SUBSTITUTION_OUTPUT_PATH}"
    )

    deletion_corrections_map = build_deletion_corrections_map(max_edits=1)
    with open(DELETION_OUTPUT_PATH, "w", encoding="utf-8") as f:
        json.dump(deletion_corrections_map, f, ensure_ascii=False, indent=2)

    print(f"\nWrote {len(deletion_corrections_map)} mappings to {DELETION_OUTPUT_PATH}")

    insertion_corrections_map = build_insertion_corrections_map(max_edits=1)
    with open(INSERTION_OUTPUT_PATH, "w", encoding="utf-8") as f:
        json.dump(insertion_corrections_map, f, ensure_ascii=False, indent=2)

    print(
        f"\nWrote {len(insertion_corrections_map)} mappings to {INSERTION_OUTPUT_PATH}"
    )

    transposition_corrections_map = build_transposition_corrections_map(max_edits=1)
    with open(TRANSPOSITION_OUTPUT_PATH, "w", encoding="utf-8") as f:
        json.dump(transposition_corrections_map, f, ensure_ascii=False, indent=2)

    print(
        f"\nWrote {len(transposition_corrections_map)} mappings to {TRANSPOSITION_OUTPUT_PATH}"
    )


if __name__ == "__main__":
    main()
