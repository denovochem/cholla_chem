import os
import sqlite3
import sys
from pathlib import Path

import pytest

# Ensure project root is on sys.path so we can import cholla_chem modules
PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

from cholla_chem.main import (  # noqa: E402
    SQLiteLookupNameResolver,
)
from cholla_chem.resolvers.lookup_resolver import (  # noqa: E402
    sqlite_lookup_resolver as slr,
)


@pytest.fixture
def sample_db(tmp_path: Path) -> Path:
    """Create a temporary SQLite database matching the lookup resolver schema."""
    db_path = tmp_path / "test_lookup.db"
    conn = sqlite3.connect(str(db_path))
    conn.executescript(
        """
        CREATE TABLE compounds (
            cid        INTEGER PRIMARY KEY,
            inchikey   TEXT,
            smiles     TEXT,
            title      TEXT,
            iupac_name TEXT
        );
        CREATE INDEX idx_compounds_inchikey ON compounds(inchikey);
        CREATE INDEX idx_compounds_smiles ON compounds(smiles);

        CREATE TABLE synonyms (
            id            INTEGER PRIMARY KEY AUTOINCREMENT,
            cid           INTEGER NOT NULL,
            synonym_text  TEXT NOT NULL,
            synonym_type  TEXT
        );
        CREATE INDEX idx_synonyms_cid ON synonyms(cid);
        CREATE INDEX idx_synonyms_text ON synonyms(synonym_text);

        CREATE VIRTUAL TABLE synonyms_fts USING fts5(
            synonym_text, content='synonyms', content_rowid='id'
        );

        CREATE TABLE md5_cid (
            md5 TEXT PRIMARY KEY,
            cid INTEGER NOT NULL
        ) WITHOUT ROWID;
        CREATE INDEX idx_md5_cid_cid ON md5_cid(cid);

        CREATE TABLE build_phases (
            phase        TEXT PRIMARY KEY,
            completed_at TEXT
        );

        INSERT INTO compounds (cid, smiles, title, iupac_name) VALUES
            (1, 'CCO', 'ethanol', 'ethanol'),
            (2, 'C', 'methane', 'methane'),
            (3, 'O', 'water', 'oxidane');

        INSERT INTO synonyms (cid, synonym_text, synonym_type) VALUES
            (1, 'ethanol', 'name'),
            (1, 'ethyl alcohol', 'synonym'),
            (2, 'methane', 'name'),
            (3, 'water', 'name'),
            (3, 'H2O', 'formula');

        INSERT INTO synonyms_fts(rowid, synonym_text)
        SELECT id, synonym_text FROM synonyms;
        """
    )
    conn.close()
    return db_path


class TestNameToSmilesSqliteLookup:
    def test_empty_list_returns_empty_dict(self, sample_db: Path):
        result = slr.name_to_smiles_sqlite_lookup([], sample_db)
        assert result == {}

    def test_exact_match_returns_smiles(self, sample_db: Path):
        result = slr.name_to_smiles_sqlite_lookup(
            ["ethanol"], sample_db, match_mode="exact"
        )
        assert result == {"ethanol": "CCO"}

    def test_exact_match_case_insensitive(self, sample_db: Path):
        result = slr.name_to_smiles_sqlite_lookup(
            ["ETHANOL"], sample_db, match_mode="exact"
        )
        assert result == {"ETHANOL": "CCO"}

    def test_fts_match_returns_smiles(self, sample_db: Path):
        result = slr.name_to_smiles_sqlite_lookup(
            ["ethanol"], sample_db, match_mode="fts"
        )
        assert result == {"ethanol": "CCO"}

    def test_no_match_returns_empty_string(self, sample_db: Path):
        result = slr.name_to_smiles_sqlite_lookup(
            ["nonexistent_compound"], sample_db, match_mode="exact"
        )
        assert result == {"nonexistent_compound": ""}

    def test_multiple_names_mixed_results(self, sample_db: Path):
        result = slr.name_to_smiles_sqlite_lookup(
            ["ethanol", "nonexistent", "water"], sample_db, match_mode="exact"
        )
        assert result["ethanol"] == "CCO"
        assert result["water"] == "O"
        assert result["nonexistent"] == ""

    def test_chunking_with_over_999_names(self, sample_db: Path):
        names = [f"compound_{i}" for i in range(1001)]
        # Insert only the first name so the rest resolve to empty strings.
        conn = sqlite3.connect(str(sample_db))
        conn.execute(
            "INSERT INTO synonyms (cid, synonym_text) VALUES (1, ?)", (names[0],)
        )
        conn.execute(
            "INSERT INTO synonyms_fts(rowid, synonym_text) SELECT id, synonym_text FROM synonyms WHERE synonym_text = ?",
            (names[0],),
        )
        conn.commit()
        conn.close()

        result = slr.name_to_smiles_sqlite_lookup(names, sample_db, match_mode="exact")
        assert result[names[0]] == "CCO"
        assert all(result[n] == "" for n in names[1:])

    def test_missing_db_returns_empty_strings(self, tmp_path: Path):
        missing_db = tmp_path / "does_not_exist.db"
        result = slr.name_to_smiles_sqlite_lookup(
            ["ethanol"], missing_db, match_mode="exact"
        )
        assert result == {"ethanol": ""}

    def test_invalid_match_mode_raises_value_error(self, sample_db: Path):
        with pytest.raises(ValueError, match="Invalid match_mode"):
            slr.name_to_smiles_sqlite_lookup(
                ["ethanol"], sample_db, match_mode="invalid"
            )

    def test_whitespace_normalization(self, sample_db: Path):
        result = slr.name_to_smiles_sqlite_lookup(
            ["  ETHANOL  "], sample_db, match_mode="exact"
        )
        assert result == {"  ETHANOL  ": "CCO"}

    def test_fts_mode_no_match(self, sample_db: Path):
        result = slr.name_to_smiles_sqlite_lookup(
            ["nonexistent"], sample_db, match_mode="fts"
        )
        assert result == {"nonexistent": ""}

    def test_duplicate_input_names(self, sample_db: Path):
        result = slr.name_to_smiles_sqlite_lookup(
            ["ethanol", "ethanol"], sample_db, match_mode="exact"
        )
        assert result == {"ethanol": "CCO"}


class TestSQLiteLookupNameResolver:
    def test_init_and_name_to_smiles(self, sample_db: Path):
        resolver = SQLiteLookupNameResolver(
            resolver_name="test_local",
            db_path=sample_db,
            resolver_weight=5.0,
            match_mode="exact",
        )
        assert resolver.resolver_name == "test_local"
        assert resolver.resolver_weight == 5.0
        assert resolver.requires_internet is False

        resolved, failures = resolver.name_to_smiles(["ethanol", "missing"])
        assert resolved["ethanol"] == "CCO"
        assert resolved["missing"] == ""
        assert failures == {}

    def test_default_match_mode_is_exact(self, sample_db: Path):
        resolver = SQLiteLookupNameResolver(
            resolver_name="test_default",
            db_path=sample_db,
        )
        assert resolver._match_mode == "exact"
