import os
import sys

# Ensure project root is on sys.path so we can import cholla_chem modules
PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

from cholla_chem.resolvers.manual_resolver import (
    _normalize_name,
    name_to_smiles_manual,
    process_name_dict,
)


def test_process_name_dict_basic_behavior():
    """process_name_dict should match names case-insensitively and strip whitespace."""
    compound_names = [" Ethanol ", "WATER", "Acetone"]
    name_dict = {
        "ethanol": "C2H6O",
        " water  ": "H2O",
        "acetone": "C3H6O",
        # Unused or empty entries should not appear in the result
        "methane": "",
    }

    result = process_name_dict(compound_names, name_dict)

    # Keys in the result should preserve the original input formatting
    assert set(result.keys()) == {" Ethanol ", "WATER", "Acetone"}
    assert result[" Ethanol "] == "C2H6O"
    assert result["WATER"] == "H2O"
    assert result["Acetone"] == "C3H6O"


def test_name_to_smiles_manual_with_provided_dict(monkeypatch):
    """name_to_smiles_manual should use the provided dict and not call the default loader."""

    # Ensure the default loader is not used when a dict is provided
    def fake_loader():  # pragma: no cover - validated by not being called
        raise AssertionError(
            "load_default_manual_name_dict should not be called when provided_name_dict is given"
        )

    monkeypatch.setattr(
        "cholla_chem.resolvers.manual_resolver.load_default_manual_name_dict",
        fake_loader,
        raising=True,
    )

    names = ["Ethanol", "Water"]
    provided_dict = {
        "ethanol": "C2H6O",
        "water": "H2O",
    }

    result = name_to_smiles_manual(names, provided_name_dict=provided_dict)

    assert result == {"Ethanol": "C2H6O", "Water": "H2O"}


def test_name_to_smiles_manual_uses_default_dict(monkeypatch):
    """When no dict is provided, name_to_smiles_manual should use load_default_manual_name_dict."""

    # Provide a small fake default dictionary via monkeypatch
    def fake_loader():
        return {
            "ethanol": "C2H6O",
            "water": "H2O",
        }

    monkeypatch.setattr(
        "cholla_chem.resolvers.manual_resolver.load_default_manual_name_dict",
        fake_loader,
        raising=True,
    )

    # Also monkeypatch logger.info so we don't depend on real logging output
    log_messages = []

    def fake_logger_info(msg):
        log_messages.append(msg)

    monkeypatch.setattr(
        "cholla_chem.resolvers.manual_resolver.logger",
        type("_L", (), {"info": staticmethod(fake_logger_info)}),
        raising=True,
    )

    names = ["Ethanol", "Water"]
    result = name_to_smiles_manual(names)

    assert result == {"Ethanol": "C2H6O", "Water": "H2O"}


def test_normalize_name():
    """_normalize_name should lowercase, strip, and remove all spaces."""
    assert _normalize_name("Na 2 SO 4") == "na2so4"
    assert _normalize_name("  CH2 Cl2  ") == "ch2cl2"
    assert _normalize_name("H2O") == "h2o"
    assert _normalize_name("  Petroleum Ether  ") == "petroleumether"


def test_process_name_dict_space_removal():
    """process_name_dict should match names with spaces to dict keys without spaces."""
    compound_names = ["Na 2 SO 4", "CH 2 Cl 2", "N 2"]
    name_dict = {
        "Na2SO4": "[Na+].[Na+].[O-]S([O-])(=O)=O",
        "CH2Cl2": "ClCCl",
        "N2": "N#N",
    }

    result = process_name_dict(compound_names, name_dict)

    assert result == {
        "Na 2 SO 4": "[Na+].[Na+].[O-]S([O-])(=O)=O",
        "CH 2 Cl 2": "ClCCl",
        "N 2": "N#N",
    }


def test_process_name_dict_duplicate_keys_prefer_nonempty():
    """When multiple dict keys normalize to the same string, non-empty values win."""
    compound_names = ["MgSO4"]
    name_dict = {
        "Mg SO4": "[Mg+2].[O-]S([O-])(=O)=O",
        "MgSO 4": "",
    }

    result = process_name_dict(compound_names, name_dict)

    assert result == {"MgSO4": "[Mg+2].[O-]S([O-])(=O)=O"}


def test_process_name_dict_duplicate_keys_empty_then_nonempty():
    """An empty value followed by a non-empty value for the same normalized key should use the non-empty one."""
    compound_names = ["MgSO4"]
    name_dict = {
        "MgSO 4": "",
        "Mg SO4": "[Mg+2].[O-]S([O-])(=O)=O",
    }

    result = process_name_dict(compound_names, name_dict)

    assert result == {"MgSO4": "[Mg+2].[O-]S([O-])(=O)=O"}


def test_process_name_dict_duplicate_keys_both_nonempty_keeps_first():
    """When two non-empty values collide, the first one encountered should be kept."""
    compound_names = ["NH4Cl"]
    name_dict = {
        "NH4 Cl": "[Cl-].[NH4+]",
        "NH 4 Cl": "NCl",
    }

    result = process_name_dict(compound_names, name_dict)

    assert result == {"NH4Cl": "[Cl-].[NH4+]"}
