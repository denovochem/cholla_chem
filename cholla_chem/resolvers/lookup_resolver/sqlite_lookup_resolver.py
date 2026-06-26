"""Resolve compound names to SMILES using a local SQLite database."""

from __future__ import annotations

import re
import sqlite3
from pathlib import Path
from typing import Dict, List, Union

from cholla_chem.utils.logging_config import logger

WS_RE = re.compile(r"\s+")
_SQLITE_MAX_VARIABLES = 999


def normalize_name(name: str) -> str:
    """
    Normalize a chemical name for consistent lookups.

    Strips leading and trailing whitespace, converts to lowercase, and collapses
    consecutive whitespace characters into a single space.

    Args:
        name (str): The raw chemical name.

    Returns:
        str: The normalized name.
    """
    name = name.strip().lower()
    name = WS_RE.sub(" ", name)
    return name


def name_to_smiles_sqlite_lookup(
    compound_name_list: List[str],
    db_path: Union[str, Path],
    match_mode: str = "exact",
) -> Dict[str, str]:
    """
    Convert chemical names to SMILES using a local SQLite database.

    Looks up each name in the database's ``synonyms`` table (joined to the
    ``compounds`` table) and returns the corresponding ``smiles`` value.
    Two matching strategies are supported:

    - ``"exact"``: Normalizes input names (lowercase, collapse whitespace) and
      queries ``LOWER(synonym_text)`` for exact matches.
    - ``"fts"``: Uses the ``synonyms_fts`` FTS5 virtual table for full-text
      synonym matching.

    Args:
        compound_name_list (List[str]): List of chemical names to resolve.
        db_path (Union[str, Path]): Path to the local SQLite database.
        match_mode (str): Matching strategy, either ``"exact"`` or ``"fts"``.
            Defaults to ``"exact"``.

    Returns:
        Dict[str, str]: Dictionary mapping each input name to its resolved SMILES
            string, or an empty string if no match was found.

    Raises:
        ValueError: If ``match_mode`` is not ``"exact"`` or ``"fts"``.
    """
    if not compound_name_list:
        return {}

    if match_mode not in ("exact", "fts"):
        raise ValueError(
            f"Invalid match_mode: {match_mode!r}. Must be 'exact' or 'fts'."
        )

    db_path = Path(db_path)
    if not db_path.exists():
        logger.warning(f"Local lookup database not found: {db_path}")
        return {name: "" for name in compound_name_list}

    try:
        conn = sqlite3.connect(str(db_path))
    except (sqlite3.Error, OSError) as e:
        logger.warning(f"Failed to open local lookup database: {e}")
        return {name: "" for name in compound_name_list}

    # Map normalized names back to their original input names.
    # Multiple originals may normalize to the same value.
    normalized_to_originals: Dict[str, List[str]] = {}
    for name in compound_name_list:
        norm = normalize_name(name)
        normalized_to_originals.setdefault(norm, []).append(name)

    results: Dict[str, str] = {}

    try:
        if match_mode == "exact":
            normalized_names = list(normalized_to_originals.keys())
            for i in range(0, len(normalized_names), _SQLITE_MAX_VARIABLES):
                chunk = normalized_names[i : i + _SQLITE_MAX_VARIABLES]
                placeholders = ",".join("?" * len(chunk))
                query = (
                    f"SELECT s.synonym_text, c.smiles "
                    f"FROM compounds c "
                    f"JOIN synonyms s ON c.cid = s.cid "
                    f"WHERE s.synonym_text IN ({placeholders})"
                )
                for db_synonym, smiles in conn.execute(query, chunk).fetchall():
                    norm_db_synonym = normalize_name(db_synonym)
                    if norm_db_synonym in normalized_to_originals:
                        for original_name in normalized_to_originals[norm_db_synonym]:
                            if original_name not in results:
                                results[original_name] = smiles
        else:  # "fts"
            for norm_name, original_names in normalized_to_originals.items():
                query = (
                    "SELECT DISTINCT c.smiles "
                    "FROM synonyms_fts fts "
                    "JOIN synonyms s ON s.id = fts.rowid "
                    "JOIN compounds c ON c.cid = s.cid "
                    "WHERE fts.synonym_text MATCH ? "
                    "LIMIT 1"
                )
                row = conn.execute(query, (norm_name,)).fetchone()
                smiles = row[0] if row else ""
                for original_name in original_names:
                    results[original_name] = smiles

    except sqlite3.Error as e:
        logger.warning(f"SQLite error during local lookup: {e}")
        return {name: "" for name in compound_name_list}
    finally:
        conn.close()

    # Fill unresolved names with empty strings so downstream logic stays stable.
    for name in compound_name_list:
        if name not in results:
            results[name] = ""

    return results
