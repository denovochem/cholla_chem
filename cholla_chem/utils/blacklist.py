"""Blacklist loading and filtering for compound names.

Provides utilities to load a list of names known to be incorrectly resolved
by database-backed resolvers (e.g. patent/paper labels like "Example 2",
"Compound 9", "IV") and filter them from compound name lists before resolution.
"""

from __future__ import annotations

import json
from importlib import resources
from typing import FrozenSet, List, Optional

from cholla_chem.utils.logging_config import logger

_BLACKLIST_CACHE: Optional[FrozenSet[str]] = None


def get_blacklist_set() -> FrozenSet[str]:
    """
    Load and return the set of blacklisted compound names as lowercase strings.

    The blacklist file is loaded once and cached for the lifetime of the process.
    All names are stored as lowercase for case-insensitive matching.

    Returns:
        FrozenSet[str]: A frozenset of lowercase blacklisted names. Returns an
        empty frozenset if the file is missing or unreadable.
    """
    global _BLACKLIST_CACHE
    if _BLACKLIST_CACHE is not None:
        return _BLACKLIST_CACHE

    try:
        raw_names = json.loads(
            resources.files("cholla_chem.datafiles")
            .joinpath("blacklisted_names.json")
            .read_text(encoding="utf-8")
        )
        _BLACKLIST_CACHE = frozenset(name.lower() for name in raw_names)
    except (FileNotFoundError, json.JSONDecodeError, OSError) as e:
        logger.warning(f"Could not load blacklisted_names.json: {e}")
        _BLACKLIST_CACHE = frozenset()

    return _BLACKLIST_CACHE


def filter_blacklisted(
    names: List[str],
    blacklist: Optional[FrozenSet[str]] = None,
) -> List[str]:
    """
    Filter out blacklisted names from a list of compound names.

    Matching is case-insensitive: each name is lowercased before checking
    membership in the blacklist set.

    Args:
        names (List[str]): List of compound names to filter.
        blacklist (Optional[FrozenSet[str]]): A pre-loaded blacklist set of
            lowercase names. If None, the default blacklist is loaded via
            `get_blacklist_set`.

    Returns:
        List[str]: A list containing only the names that are not blacklisted,
        preserving the original casing and order of the input.
    """
    if blacklist is None:
        blacklist = get_blacklist_set()
    if not blacklist:
        return list(names)
    return [name for name in names if name.lower() not in blacklist]
