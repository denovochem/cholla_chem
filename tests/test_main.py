import os
import sys
from typing import Dict, List, Tuple

import pytest

PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

from cholla_chem.main import (
    ChemicalNameResolver,
    assemble_compounds_resolution_dict,
    assemble_split_compounds_resolution_dict,
    get_resolvers_weight_dict,
    resolve_compounds_to_smiles,
    resolve_compounds_using_resolvers,
    select_smiles_with_criteria,
)


class DummyResolver(ChemicalNameResolver):
    """In-memory resolver for testing. Tracks which compounds it was called with."""

    def __init__(
        self,
        resolver_name: str,
        mapping: Dict[str, str],
        info_mapping: Dict[str, str] | None = None,
        weight: float = 1.0,
    ):
        super().__init__("dummy", resolver_name, weight)
        self._mapping = mapping
        self._info_mapping = info_mapping or {}
        self.call_log: List[List[str]] = []

    def name_to_smiles(
        self, compound_name_list: List[str]
    ) -> Tuple[Dict[str, str], Dict[str, str]]:
        self.call_log.append(list(compound_name_list))
        out = {name: self._mapping.get(name, "") for name in compound_name_list}
        info = {name: self._info_mapping.get(name, "") for name in compound_name_list}
        info = {k: v for k, v in info.items() if v}
        return out, info


# ---------------------------------------------------------------------------
# Baseline tests (restore coverage of existing behavior)
# ---------------------------------------------------------------------------


def test_get_resolvers_weight_dict_uses_resolver_name_and_weight():
    r1 = DummyResolver("r1", {}, weight=1.5)
    r2 = DummyResolver("r2", {}, weight=3.0)

    result = get_resolvers_weight_dict([r1, r2])

    assert result == {"r1": 1.5, "r2": 3.0}


def test_resolve_compounds_using_resolvers_basic():
    r1 = DummyResolver("res1", {"a": "S_a1", "b": "S_b1"}, {"a": "i1"})
    r2 = DummyResolver("res2", {"a": "S_a2"}, {"b": "i2"})

    out = resolve_compounds_using_resolvers(["a", "b"], [r1, r2], batch_size=500)

    assert set(out.keys()) == {"res1", "res2"}
    assert out["res1"]["out"] == {"a": "S_a1", "b": "S_b1"}
    assert out["res1"]["additional_info"] == {"a": "i1"}
    assert out["res2"]["out"] == {"a": "S_a2", "b": ""}
    assert out["res2"]["additional_info"] == {"b": "i2"}


def test_resolve_compounds_using_resolvers_batches():
    r1 = DummyResolver("res1", {"a": "S_a", "b": "S_b"})

    out = resolve_compounds_using_resolvers(["a", "b"], [r1], batch_size=1)

    assert out["res1"]["out"] == {"a": "S_a", "b": "S_b"}
    assert len(r1.call_log) == 2
    assert r1.call_log[0] == ["a"]
    assert r1.call_log[1] == ["b"]


def test_assemble_compounds_resolution_dict_uses_canonical_smiles(monkeypatch):
    def fake_canonicalize(smiles: str) -> str:
        return smiles.upper() if smiles else ""

    monkeypatch.setattr(
        "cholla_chem.main.canonicalize_smiles",
        fake_canonicalize,
        raising=True,
    )

    compounds = ["ethanol"]
    cleaned_mapping = {"ethanol": "ethanol_clean"}
    resolvers_out = {
        "r1": {"out": {"ethanol_clean": "c2h6o"}, "additional_info": {}},
        "r2": {
            "out": {"ethanol_clean": "C2H6O"},
            "additional_info": {"ethanol_clean": "ok"},
        },
    }

    result = assemble_compounds_resolution_dict(
        compounds, resolvers_out, cleaned_mapping
    )

    assert set(result.keys()) == {"ethanol"}
    entry = result["ethanol"]
    assert entry["SMILES"] == ""
    assert entry["SMILES_dict"] == {"C2H6O": ["r1", "r2"]}
    assert entry["additional_info"] == {"r2": "ok"}


def test_assemble_split_compounds_resolution_dict_merges_split_smiles(monkeypatch):
    monkeypatch.setattr(
        "cholla_chem.main.canonicalize_smiles",
        lambda s: s,
        raising=True,
    )

    resolvers_out_dict = {}
    cleaned_mapping = {"mix": "mix_clean"}
    delimiter_split_dict = {"mix_clean": ["a", "b"]}
    compounds_out_dict = {
        "mix": {
            "SMILES": "",
            "SMILES_source": [],
            "SMILES_dict": {},
            "additional_info": {},
        }
    }

    def fake_resolve_delimiter_split_dict(compound_cleaned, ro_dict, del_dict):
        assert compound_cleaned == "mix_clean"
        assert del_dict is delimiter_split_dict
        return {
            "S_mix1": ["r_split1"],
            "S_mix2": ["r_split2", "r_split3"],
        }

    monkeypatch.setattr(
        "cholla_chem.main.resolve_delimiter_split_dict",
        fake_resolve_delimiter_split_dict,
        raising=True,
    )

    result = assemble_split_compounds_resolution_dict(
        compounds_out_dict,
        ["mix"],
        resolvers_out_dict,
        cleaned_mapping,
        delimiter_split_dict,
    )

    entry = result["mix"]
    assert entry["SMILES_dict"] == {
        "S_mix1": ["r_split1"],
        "S_mix2": ["r_split2", "r_split3"],
    }


def test_select_smiles_with_criteria_uses_selector(monkeypatch):
    compounds_out = {
        "ethanol": {
            "SMILES": "",
            "SMILES_source": [],
            "SMILES_dict": {
                "S1": ["r1"],
                "S2": ["r2"],
            },
            "additional_info": {},
        }
    }

    resolvers_weight = {"r1": 1.0, "r2": 2.0}
    priority_order = ["r2", "r1"]

    selected_calls: List[Tuple] = []

    class FakeSelector:
        def __init__(self, compounds_out_dict, weight_dict, priority_list):
            selected_calls.append((compounds_out_dict, weight_dict, priority_list))

        def select_smiles(self, compound, mode):
            assert compound == "ethanol"
            assert mode == "weighted"
            return "S2", ["r2"]

    monkeypatch.setattr(
        "cholla_chem.main.SMILESSelector",
        FakeSelector,
        raising=True,
    )

    result = select_smiles_with_criteria(
        compounds_out,
        resolvers_weight,
        priority_order,
        smiles_selection_mode="weighted",
    )

    assert selected_calls
    assert result["ethanol"]["SMILES"] == "S2"
    assert result["ethanol"]["SMILES_source"] == ["r2"]


def test_resolve_compounds_to_smiles_validates_inputs():
    with pytest.raises(ValueError):
        resolve_compounds_to_smiles(123, resolvers_list=[DummyResolver("r", {})])

    with pytest.raises(ValueError):
        resolve_compounds_to_smiles([], resolvers_list=[DummyResolver("r", {})])

    with pytest.raises(ValueError):
        resolve_compounds_to_smiles(["ok", 1], resolvers_list=[DummyResolver("r", {})])

    class NotAResolver:
        pass

    with pytest.raises(ValueError):
        resolve_compounds_to_smiles(["a"], resolvers_list=[NotAResolver()])

    r1 = DummyResolver("dup", {})
    r2 = DummyResolver("dup", {})
    with pytest.raises(ValueError):
        resolve_compounds_to_smiles(["a"], resolvers_list=[r1, r2])

    with pytest.raises(ValueError):
        resolve_compounds_to_smiles(
            ["a"],
            resolvers_list=[DummyResolver("r", {})],
            smiles_selection_mode=123,
        )

    with pytest.raises(ValueError):
        resolve_compounds_to_smiles(
            ["a"],
            resolvers_list=[DummyResolver("r", {})],
            detailed_name_dict="yes",
        )

    with pytest.raises(TypeError):
        resolve_compounds_to_smiles(
            ["a"],
            resolvers_list=[DummyResolver("r", {})],
            batch_size=1.5,
        )

    with pytest.raises(ValueError):
        resolve_compounds_to_smiles(
            ["a"],
            resolvers_list=[DummyResolver("r", {})],
            batch_size=0,
        )

    with pytest.raises(ValueError):
        resolve_compounds_to_smiles(
            ["a"],
            resolvers_list=[DummyResolver("r", {})],
            batch_size=1001,
        )

    with pytest.raises(ValueError):
        resolve_compounds_to_smiles(
            ["a"],
            resolvers_list=[DummyResolver("r", {})],
            split_names_to_solve="yes",
        )


def test_resolve_compounds_to_smiles_happy_path(monkeypatch):
    resolver = DummyResolver("dummy_res", {"a": "S_a"})

    monkeypatch.setattr(
        "cholla_chem.main.normalize_unicode_and_return_mapping",
        lambda names: (names, {n: n for n in names}),
        raising=True,
    )
    monkeypatch.setattr(
        "cholla_chem.main.split_compounds_on_delimiters_and_return_mapping",
        lambda names: (names, {}),
        raising=True,
    )
    monkeypatch.setattr(
        "cholla_chem.utils.chem_utils.canonicalize_smiles",
        lambda s: s,
        raising=True,
    )

    class FakeSelector:
        def __init__(self, compounds_out_dict, weight_dict, priority_list):
            self._d = compounds_out_dict

        def select_smiles(self, compound, mode):
            entry = self._d[compound]
            smiles = next(iter(entry["SMILES_dict"].keys()), "")
            return smiles, entry["SMILES_dict"].get(smiles, [])

    monkeypatch.setattr(
        "cholla_chem.main.SMILESSelector",
        FakeSelector,
        raising=True,
    )

    result = resolve_compounds_to_smiles(["a"], resolvers_list=[resolver])

    assert result == {"a": "S_a"}


# ---------------------------------------------------------------------------
# Early exit tests
# ---------------------------------------------------------------------------


def test_exit_early_skips_resolved_compounds(monkeypatch):
    monkeypatch.setattr(
        "cholla_chem.main.canonicalize_smiles",
        lambda s: s,
        raising=True,
    )

    r1 = DummyResolver("r1", {"a": "S_a"})
    r2 = DummyResolver("r2", {"b": "S_b"})

    out = resolve_compounds_using_resolvers(
        ["a", "b"], [r1, r2], batch_size=500, exit_early=True
    )

    assert out["r1"]["out"] == {"a": "S_a", "b": ""}
    assert "a" not in out["r2"]["out"]
    assert out["r2"]["out"] == {"b": "S_b"}
    assert r2.call_log == [["b"]]


def test_exit_early_all_resolved_breaks_loop(monkeypatch):
    monkeypatch.setattr(
        "cholla_chem.main.canonicalize_smiles",
        lambda s: s,
        raising=True,
    )

    r1 = DummyResolver("r1", {"a": "S_a", "b": "S_b"})
    r2 = DummyResolver("r2", {})

    out = resolve_compounds_using_resolvers(
        ["a", "b"], [r1, r2], batch_size=500, exit_early=True
    )

    assert out["r1"]["out"] == {"a": "S_a", "b": "S_b"}
    assert r2.call_log == []


def test_exit_early_none_resolved_calls_all_resolvers(monkeypatch):
    monkeypatch.setattr(
        "cholla_chem.main.canonicalize_smiles",
        lambda s: s,
        raising=True,
    )

    r1 = DummyResolver("r1", {})
    r2 = DummyResolver("r2", {"a": "S_a"})

    out = resolve_compounds_using_resolvers(
        ["a"], [r1, r2], batch_size=500, exit_early=True
    )

    assert out["r1"]["out"] == {"a": ""}
    assert out["r2"]["out"] == {"a": "S_a"}
    assert r1.call_log == [["a"]]
    assert r2.call_log == [["a"]]


def test_exit_early_invalid_smiles_not_counted_as_resolved(monkeypatch):
    def fake_canonicalize(smiles: str) -> str:
        if smiles == "invalid":
            return ""
        return smiles

    monkeypatch.setattr(
        "cholla_chem.main.canonicalize_smiles",
        fake_canonicalize,
        raising=True,
    )

    r1 = DummyResolver("r1", {"a": "invalid"})
    r2 = DummyResolver("r2", {"a": "S_a"})

    out = resolve_compounds_using_resolvers(
        ["a"], [r1, r2], batch_size=500, exit_early=True
    )

    assert out["r1"]["out"] == {"a": "invalid"}
    assert out["r2"]["out"] == {"a": "S_a"}
    assert r2.call_log == [["a"]]


def test_exit_early_preserves_output_structure(monkeypatch):
    monkeypatch.setattr(
        "cholla_chem.main.canonicalize_smiles",
        lambda s: s,
        raising=True,
    )

    r1 = DummyResolver("r1", {"a": "S_a"}, {"a": "info_a"})
    r2 = DummyResolver("r2", {"b": "S_b"})

    out = resolve_compounds_using_resolvers(
        ["a", "b"], [r1, r2], batch_size=500, exit_early=True
    )

    assert set(out.keys()) == {"r1", "r2"}
    assert set(out["r1"].keys()) == {"out", "additional_info"}
    assert set(out["r2"].keys()) == {"out", "additional_info"}
    assert out["r1"]["additional_info"] == {"a": "info_a"}
    assert out["r2"]["additional_info"] == {}


def test_exit_early_end_to_end(monkeypatch):
    monkeypatch.setattr(
        "cholla_chem.main.normalize_unicode_and_return_mapping",
        lambda names: (names, {n: n for n in names}),
        raising=True,
    )
    monkeypatch.setattr(
        "cholla_chem.main.split_compounds_on_delimiters_and_return_mapping",
        lambda names: (names, {}),
        raising=True,
    )
    monkeypatch.setattr(
        "cholla_chem.main.canonicalize_smiles",
        lambda s: s,
        raising=True,
    )

    class FakeSelector:
        def __init__(self, compounds_out_dict, weight_dict, priority_list):
            self._d = compounds_out_dict

        def select_smiles(self, compound, mode):
            entry = self._d[compound]
            smiles = next(iter(entry["SMILES_dict"].keys()), "")
            return smiles, entry["SMILES_dict"].get(smiles, [])

    monkeypatch.setattr(
        "cholla_chem.main.SMILESSelector",
        FakeSelector,
        raising=True,
    )

    r1 = DummyResolver("r1", {"aspirin": "CC(=O)Oc1ccccc1C(=O)O"})
    r2 = DummyResolver("r2", {"aspirin": "OTHER_SMILES"})
    r3 = DummyResolver("r3", {})

    result = resolve_compounds_to_smiles(
        ["aspirin"], resolvers_list=[r1, r2, r3], exit_early=True
    )

    assert result == {"aspirin": "CC(=O)Oc1ccccc1C(=O)O"}
    assert r1.call_log == [["aspirin"]]
    assert r2.call_log == []
    assert r3.call_log == []


def test_exit_early_validation():
    with pytest.raises(ValueError):
        resolve_compounds_to_smiles(
            ["a"],
            resolvers_list=[DummyResolver("r", {})],
            exit_early="yes",
        )


def test_exit_early_passes_to_recursive_name_correction(monkeypatch):
    monkeypatch.setattr(
        "cholla_chem.main.normalize_unicode_and_return_mapping",
        lambda names: (names, {n: n for n in names}),
        raising=True,
    )
    monkeypatch.setattr(
        "cholla_chem.main.split_compounds_on_delimiters_and_return_mapping",
        lambda names: (names, {}),
        raising=True,
    )
    monkeypatch.setattr(
        "cholla_chem.main.canonicalize_smiles",
        lambda s: s,
        raising=True,
    )

    class FakeSelector:
        def __init__(self, compounds_out_dict, weight_dict, priority_list):
            self._d = compounds_out_dict

        def select_smiles(self, compound, mode):
            entry = self._d[compound]
            smiles = next(iter(entry["SMILES_dict"].keys()), "")
            return smiles, entry["SMILES_dict"].get(smiles, [])

    monkeypatch.setattr(
        "cholla_chem.main.SMILESSelector",
        FakeSelector,
        raising=True,
    )

    r1 = DummyResolver("r1", {"corrected_name": "S_corrected"})

    resolve_calls: List[bool] = []

    original_resolve = resolve_compounds_using_resolvers

    def spy_resolve(compounds_list, resolvers_list, batch_size, exit_early=False):
        resolve_calls.append(exit_early)
        return original_resolve(compounds_list, resolvers_list, batch_size, exit_early)

    monkeypatch.setattr(
        "cholla_chem.main.resolve_compounds_using_resolvers",
        spy_resolve,
        raising=True,
    )

    monkeypatch.setattr(
        "cholla_chem.main.correct_names",
        lambda compounds_out_dict, config, resolve_peptide: {
            "unresolved_compound": {
                "selected_name": "corrected_name",
                "top_5": [],
                "name_manipulation_method": "test",
                "SMILES": "",
            }
        },
        raising=True,
    )

    result = resolve_compounds_to_smiles(
        ["unresolved_compound"],
        resolvers_list=[r1],
        exit_early=True,
        detailed_name_dict=True,
    )

    assert len(resolve_calls) == 2
    assert all(call is True for call in resolve_calls)
    assert result["unresolved_compound"]["SMILES"] == "S_corrected"


# ---------------------------------------------------------------------------
# Blacklist integration tests
# ---------------------------------------------------------------------------


def _setup_blacklist_mocks(monkeypatch):
    """Apply the standard mocks needed for resolve_compounds_to_smiles tests."""
    monkeypatch.setattr(
        "cholla_chem.main.normalize_unicode_and_return_mapping",
        lambda names: (names, {n: n for n in names}),
        raising=True,
    )
    monkeypatch.setattr(
        "cholla_chem.main.split_compounds_on_delimiters_and_return_mapping",
        lambda names: (names, {}),
        raising=True,
    )
    monkeypatch.setattr(
        "cholla_chem.utils.chem_utils.canonicalize_smiles",
        lambda s: s,
        raising=True,
    )

    class FakeSelector:
        def __init__(self, compounds_out_dict, weight_dict, priority_list):
            self._d = compounds_out_dict

        def select_smiles(self, compound, mode):
            entry = self._d[compound]
            smiles = next(iter(entry["SMILES_dict"].keys()), "")
            return smiles, entry["SMILES_dict"].get(smiles, [])

    monkeypatch.setattr(
        "cholla_chem.main.SMILESSelector",
        FakeSelector,
        raising=True,
    )


def test_blacklisted_name_gets_empty_smiles(monkeypatch):
    """A blacklisted name should appear in the output with empty SMILES."""
    _setup_blacklist_mocks(monkeypatch)
    resolver = DummyResolver("dummy_res", {"ethanol": "S_ethanol"})

    result = resolve_compounds_to_smiles(
        ["ethanol", "Example 2"], resolvers_list=[resolver]
    )

    assert result["ethanol"] == "S_ethanol"
    assert result["Example 2"] == ""


def test_blacklisted_name_not_sent_to_resolvers(monkeypatch):
    """Blacklisted names should not be passed to any resolver."""
    _setup_blacklist_mocks(monkeypatch)
    resolver = DummyResolver("dummy_res", {"ethanol": "S_ethanol", "Example 2": "S_bad"})

    result = resolve_compounds_to_smiles(
        ["ethanol", "Example 2"], resolvers_list=[resolver]
    )

    # The resolver should only have been called with non-blacklisted names
    assert len(resolver.call_log) == 1
    assert "Example 2" not in resolver.call_log[0]
    assert "ethanol" in resolver.call_log[0]
    # Even though the resolver has a mapping for "Example 2", it should not be used
    assert result["Example 2"] == ""


def test_blacklist_case_insensitive(monkeypatch):
    """Uppercase variants of blacklisted names should also be filtered."""
    _setup_blacklist_mocks(monkeypatch)
    resolver = DummyResolver("dummy_res", {"ethanol": "S_ethanol", "COMPOUND 9": "S_bad"})

    result = resolve_compounds_to_smiles(
        ["ethanol", "COMPOUND 9"], resolvers_list=[resolver]
    )

    assert result["ethanol"] == "S_ethanol"
    assert result["COMPOUND 9"] == ""


def test_use_blacklist_false_disables_filtering(monkeypatch):
    """With use_blacklist=False, blacklisted names should be sent to resolvers."""
    _setup_blacklist_mocks(monkeypatch)
    resolver = DummyResolver("dummy_res", {"ethanol": "S_ethanol", "Example 2": "S_bad"})

    result = resolve_compounds_to_smiles(
        ["ethanol", "Example 2"], resolvers_list=[resolver], use_blacklist=False
    )

    assert result["ethanol"] == "S_ethanol"
    assert result["Example 2"] == "S_bad"


def test_blacklisted_name_excluded_from_name_correction_results(monkeypatch):
    """Even if correct_names returns a correction for a blacklisted name, it should not be applied."""
    _setup_blacklist_mocks(monkeypatch)
    resolver = DummyResolver("dummy_res", {"ethanol": "S_ethanol"})

    def fake_correct_names(compounds_out_dict, config, resolve_peptide):
        return {
            "Example 2": {
                "selected_name": "corrected_example2",
                "top_5": [],
                "name_manipulation_method": "test",
                "SMILES": "",
            },
            "ethanol": {
                "selected_name": "corrected_ethanol",
                "top_5": [],
                "name_manipulation_method": "test",
                "SMILES": "",
            },
        }

    monkeypatch.setattr(
        "cholla_chem.main.correct_names",
        fake_correct_names,
        raising=True,
    )

    # Add a resolver for the corrected names
    resolver2 = DummyResolver(
        "dummy_res2",
        {"corrected_ethanol": "S_corrected_ethanol", "corrected_example2": "S_bad"},
    )

    result = resolve_compounds_to_smiles(
        ["ethanol", "Example 2"],
        resolvers_list=[resolver, resolver2],
        detailed_name_dict=True,
    )

    # "ethanol" should get corrected and resolved
    assert result["ethanol"]["SMILES"] == "S_corrected_ethanol"
    # "Example 2" should NOT get corrected - should remain empty
    assert result["Example 2"]["SMILES"] == ""
