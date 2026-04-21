import os
import sys
from urllib.error import HTTPError, URLError

import pytest

# Ensure project root is on sys.path so we can import cholla_chem modules
PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

from cholla_chem.resolvers.pubchem_resolver import pubchem_resolver as pcr  # noqa: E402


class _FakeHTTPResponse:
    def __init__(self, body: bytes):
        self._body = body

    def read(self) -> bytes:
        return self._body

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        return False


def test_request_smiles_by_name_builds_expected_url(monkeypatch):
    captured = {}

    def fake_urlopen(api_url, timeout=None, context=None):
        captured["api_url"] = api_url
        captured["timeout"] = timeout
        return _FakeHTTPResponse(b'{"ok": true}')

    monkeypatch.setattr(pcr, "urlopen", fake_urlopen, raising=True)

    response = pcr.request_smiles_by_name("sodium chloride", timeout=7)

    assert response == b'{"ok": true}'
    assert (
        captured["api_url"]
        == f"{pcr.API_BASE}/compound/name/sodium%20chloride/property/SMILES/JSON"
    )
    assert captured["timeout"] == 7


def test_request_smiles_by_name_raises_for_empty_name():
    with pytest.raises(ValueError, match="compound_name cannot be empty"):
        pcr.request_smiles_by_name("")


def test_request_smiles_by_name_wraps_http_error(monkeypatch):
    def fake_urlopen(api_url, timeout=None, context=None):
        raise HTTPError(api_url, 404, "Not Found", hdrs=None, fp=None)

    monkeypatch.setattr(pcr, "urlopen", fake_urlopen, raising=True)

    with pytest.raises(Exception, match="HTTP Error 404: Not Found"):
        pcr.request_smiles_by_name("ethanol")


def test_request_smiles_by_name_wraps_url_error(monkeypatch):
    def fake_urlopen(api_url, timeout=None, context=None):
        raise URLError("network down")

    monkeypatch.setattr(pcr, "urlopen", fake_urlopen, raising=True)

    with pytest.raises(Exception, match="URL Error: network down"):
        pcr.request_smiles_by_name("ethanol")


def test_get_smiles_by_name_returns_smiles_from_payload(monkeypatch):
    monkeypatch.setattr(
        pcr,
        "request_smiles_by_name",
        lambda compound_name: b'{"PropertyTable": {"Properties": [{"SMILES": "CCO"}]}}',
        raising=True,
    )

    assert pcr.get_smiles_by_name("ethanol") == "CCO"


def test_get_smiles_by_name_returns_empty_string_when_properties_missing(monkeypatch):
    monkeypatch.setattr(
        pcr,
        "request_smiles_by_name",
        lambda compound_name: b'{"PropertyTable": {"Properties": []}}',
        raising=True,
    )

    assert pcr.get_smiles_by_name("ethanol") == ""


def test_get_smiles_by_name_returns_empty_string_on_request_error(monkeypatch):
    monkeypatch.setattr(
        pcr,
        "request_smiles_by_name",
        lambda compound_name: (_ for _ in ()).throw(RuntimeError("boom")),
        raising=True,
    )

    assert pcr.get_smiles_by_name("ethanol") == ""


def test_name_to_smiles_pubchem_maps_each_filtered_name(monkeypatch):
    def fake_filter_latin1_compatible(names):
        return names

    monkeypatch.setattr(
        pcr,
        "filter_latin1_compatible",
        fake_filter_latin1_compatible,
        raising=True,
    )
    monkeypatch.setattr(
        pcr,
        "get_smiles_by_name",
        lambda compound_name: f"SMILES_{compound_name}",
        raising=True,
    )

    compound_names = ["ethanol", "water", "acetone"]
    result = pcr.name_to_smiles_pubchem(compound_names)

    assert result == {
        "ethanol": "SMILES_ethanol",
        "water": "SMILES_water",
        "acetone": "SMILES_acetone",
    }


def test_name_to_smiles_pubchem_returns_empty_dict_when_filtered_names_empty(
    monkeypatch,
):
    monkeypatch.setattr(pcr, "filter_latin1_compatible", lambda names: [], raising=True)

    assert pcr.name_to_smiles_pubchem(["ethanol"]) == {}
