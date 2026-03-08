from __future__ import annotations

import json
from pathlib import Path
from typing import Callable, Dict

import pytest

from cholla_chem.name_manipulation.name_correction.build_flashtext_ocr_map import (
    build_deletion_corrections_map,
    build_insertion_corrections_map,
    build_substitution_corrections_map,
    build_transposition_corrections_map,
)

# Paths
PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATAFILES_DIR = PROJECT_ROOT / "cholla_chem" / "datafiles" / "flashtext_jsons"

SUBSTITUTIONS_MAP_PATH = DATAFILES_DIR / "substitutions_map.json"
DELETIONS_MAP_PATH = DATAFILES_DIR / "deletions_map.json"
INSERTIONS_MAP_PATH = DATAFILES_DIR / "insertions_map.json"
TRANSPOSITIONS_MAP_PATH = DATAFILES_DIR / "transpositions_map.json"


MapBuilder = Callable[[int], Dict[str, str]]


def _load_json(path: Path) -> Dict[str, str]:
    with open(path, "r", encoding="utf-8") as f:
        data = json.load(f)
    if not isinstance(data, dict):
        raise TypeError(f"Expected dict in {path}, got {type(data)}")
    return data


def _assert_map_matches_generated(
    shipped_path: Path, builder: MapBuilder, max_edits: int = 1
) -> None:
    if not shipped_path.exists():
        pytest.skip(f"{shipped_path.name} does not exist")

    shipped_map = _load_json(shipped_path)
    generated_map = builder(max_edits)

    shipped_keys = set(shipped_map.keys())
    generated_keys = set(generated_map.keys())

    missing_in_shipped = generated_keys - shipped_keys
    extra_in_shipped = shipped_keys - generated_keys

    assert not missing_in_shipped, (
        f"{shipped_path.name} is missing {len(missing_in_shipped)} keys. "
        "Run: python scripts/build_ocr_map.py"
    )

    assert not extra_in_shipped, (
        f"{shipped_path.name} has {len(extra_in_shipped)} extra keys. "
        "Run: python scripts/build_ocr_map.py"
    )

    mismatches = [k for k in shipped_keys if shipped_map[k] != generated_map[k]]
    assert not mismatches, (
        f"{shipped_path.name} has {len(mismatches)} value mismatches. "
        "Run: python scripts/build_ocr_map.py"
    )


@pytest.mark.parametrize(
    ("shipped_path", "builder"),
    [
        (SUBSTITUTIONS_MAP_PATH, build_substitution_corrections_map),
        (DELETIONS_MAP_PATH, build_deletion_corrections_map),
        (INSERTIONS_MAP_PATH, build_insertion_corrections_map),
        (TRANSPOSITIONS_MAP_PATH, build_transposition_corrections_map),
    ],
)
def test_flashtext_maps_exist_and_match_generated(
    shipped_path: Path, builder: MapBuilder
) -> None:
    assert shipped_path.exists(), (
        f"{shipped_path.name} not found at {shipped_path}. "
        "Run: python scripts/build_ocr_map.py"
    )

    _assert_map_matches_generated(shipped_path, builder, max_edits=1)
