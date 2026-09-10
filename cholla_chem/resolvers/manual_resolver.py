import json
from importlib import resources
from typing import Dict, List

from cholla_chem.utils.logging_config import logger


def load_default_manual_name_dict() -> Dict[str, str]:
    """Load manual name dictionary from package data using importlib.resources."""
    # Open the file as text from within the package
    with resources.open_text(
        "cholla_chem.datafiles.name_dicts", "manual_name_dict.json"
    ) as f:
        return json.load(f)


def _normalize_name(name: str) -> str:
    """
    Normalize a compound name for case-sensitive matching by stripping whitespace and removing spaces.

    Removing spaces handles patent-derived names where subscripts are separated by spaces
    (e.g. ``"Na 2 SO 4"`` normalizes to ``"Na2SO4"`` to match the dict key ``"Na2SO4"``).

    Unlike :func:`_normalize_name_ci`, this function preserves character casing so that
    chemically distinct compounds that differ only in case (e.g. ``"CoCl2`` vs
    ``"COCl2"``) remain distinguishable.

    Args:
        name (str): The compound name to normalize.

    Returns:
        str: The normalized name with leading/trailing whitespace removed and all
        internal spaces removed.
    """
    return name.strip().replace(" ", "")


def _normalize_name_ci(name: str) -> str:
    """
    Normalize a compound name for case-insensitive matching.

    Applies the same normalization as :func:`_normalize_name` (strip + remove spaces)
    and additionally lowercases the result. This is used for the case-insensitive
    fallback lookup when an exact case match is not found.

    Args:
        name (str): The compound name to normalize.

    Returns:
        str: The normalized, lowercased name.
    """
    return name.lower().strip().replace(" ", "")


def process_name_dict(
    compound_name_list: List[str],
    name_dict: Dict[str, str],
) -> Dict[str, str]:
    """
    Match compound names to SMILES using a two-phase lookup with case disambiguation.

    Each input name is first matched case-sensitively (after stripping whitespace and
    removing internal spaces). If no exact match is found, a case-insensitive fallback
    is attempted. When the case-insensitive key maps to exactly one case variant, the
    match is accepted. When multiple case variants exist with different SMILES, the
    name is considered ambiguous and skipped (a warning is logged).

    When multiple dict keys normalize to the same case-preserving string, non-empty
    SMILES values are preferred over empty ones, and the first non-empty value is kept.

    Args:
        compound_name_list (List[str]): A list of compound names to process.
        name_dict (Dict[str, str]): A dictionary of compound names to SMILES strings.

    Returns:
        Dict[str, str]: A dictionary mapping original input compound names to their
        resolved SMILES strings. Names that are ambiguous or not found are omitted.
    """
    # Build exact-match index: case-preserving normalized key -> SMILES.
    # When multiple dict keys collapse to the same case-preserving key, prefer
    # the first non-empty SMILES value (same behavior as the previous implementation).
    exact_index: Dict[str, str] = {}
    for k, v in name_dict.items():
        key = _normalize_name(k)
        if key not in exact_index or (not exact_index[key] and v):
            exact_index[key] = v

    # Build case-insensitive index: lowercased normalized key -> list of unique
    # case-preserving keys (only for entries with non-empty SMILES).
    ci_index: Dict[str, List[str]] = {}
    for k, v in name_dict.items():
        if not v:
            continue
        ci_key = _normalize_name_ci(k)
        cp_key = _normalize_name(k)
        if ci_key not in ci_index:
            ci_index[ci_key] = []
        if cp_key not in ci_index[ci_key]:
            ci_index[ci_key].append(cp_key)

    processed_name_dict: Dict[str, str] = {}
    for ele in compound_name_list:
        # Phase 1: exact case match.
        exact_key = _normalize_name(ele)
        if exact_index.get(exact_key):
            processed_name_dict[ele] = exact_index[exact_key]
            continue

        # Phase 2: case-insensitive fallback with disambiguation.
        ci_key = _normalize_name_ci(ele)
        if ci_key in ci_index:
            candidates = ci_index[ci_key]
            if len(candidates) == 1:
                # Unambiguous — accept the single case variant.
                processed_name_dict[ele] = exact_index[candidates[0]]
            else:
                # Ambiguous — multiple case variants with potentially different SMILES.
                logger.warning(
                    "Ambiguous compound name %r: case-insensitive match resolves to "
                    "multiple entries with different casing: %s. Skipping.",
                    ele,
                    candidates,
                )

    return processed_name_dict


def name_to_smiles_manual(
    compound_name_list: List[str],
    provided_name_dict: Dict[str, str] | None = None,
) -> Dict[str, str]:
    """
    Convert a list of compound names to their corresponding SMILES strings using a manual name dictionary.

    Args:
        compound_name_list (List[str]): A list of compound names to convert.
        provided_name_dict (Dict[str, str], optional): A manual name dictionary to use. Defaults to None.

    Returns:
        Dict[str, str]: A dictionary of converted compound names to SMILES strings.
    """
    manual_name_dict = {}
    if not provided_name_dict:
        loaded_manual_name_dict = load_default_manual_name_dict()
    else:
        loaded_manual_name_dict = provided_name_dict
        logger.info("Using provided name dictionary.")

    manual_name_dict = process_name_dict(compound_name_list, loaded_manual_name_dict)

    return manual_name_dict
