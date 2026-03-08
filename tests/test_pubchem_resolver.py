import json
import os
import sys

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


def test_name_to_smiles_pubchem_basic_mapping(monkeypatch):
    """name_to_smiles_pubchem should map each name to the SMILES returned by PubChem."""

    # Make filter_latin1_compatible a pass-through so we control the argument into get_compounds
    def fake_filter_latin1_compatible(names):
        return names

    monkeypatch.setattr(
        "cholla_chem.resolvers.pubchem_resolver.pubchem_resolver.filter_latin1_compatible",
        fake_filter_latin1_compatible,
        raising=True,
    )

    captured_args = {}

    def fake_get_compounds(names, identifier_type):
        captured_args["names"] = names
        captured_args["identifier_type"] = identifier_type
        return [f"SMILES_{name}" for name in names]

    monkeypatch.setattr(
        "cholla_chem.resolvers.pubchem_resolver.pubchem_resolver.get_compounds",
        fake_get_compounds,
        raising=True,
    )

    compound_names = ["ethanol", "water", "acetone"]
    result = pcr.name_to_smiles_pubchem(compound_names)

    # Ensure get_compounds was called with the expected arguments
    assert captured_args["names"] == compound_names
    assert captured_args["identifier_type"] == "name"

    # Ensure the result dictionary contains all names with the expected SMILES
    assert set(result.keys()) == set(compound_names)
    for name in compound_names:
        assert result[name] == f"SMILES_{name}"


def test_name_to_smiles_pubchem_handles_none_results(monkeypatch):
    """If PubChem returns None for an entry, the corresponding SMILES should be an empty string."""

    def fake_filter_latin1_compatible(names):
        return names

    monkeypatch.setattr(
        "cholla_chem.resolvers.pubchem_resolver.pubchem_resolver.filter_latin1_compatible",
        fake_filter_latin1_compatible,
        raising=True,
    )

    def fake_get_compounds(names, identifier_type):
        # Return a list with one valid compound and one None
        return ["SMILES_valid", None]

    monkeypatch.setattr(
        "cholla_chem.resolvers.pubchem_resolver.pubchem_resolver.get_compounds",
        fake_get_compounds,
        raising=True,
    )

    compound_names = ["valid_name", "missing_name"]
    result = pcr.name_to_smiles_pubchem(compound_names)

    assert result["valid_name"] == "SMILES_valid"
    assert result["missing_name"] == ""


def test_name_to_smiles_pubchem_logs_length_mismatch(monkeypatch):
    """If lengths mismatch, a warning should be logged and an empty dictionary should be returned."""

    def fake_filter_latin1_compatible(names):
        return names

    monkeypatch.setattr(
        "cholla_chem.resolvers.pubchem_resolver.pubchem_resolver.filter_latin1_compatible",
        fake_filter_latin1_compatible,
        raising=True,
    )

    # Return fewer results than input names to trigger the warning
    def fake_get_compounds(names, identifier_type):
        return ["SMILES_only_first"]

    monkeypatch.setattr(
        "cholla_chem.resolvers.pubchem_resolver.pubchem_resolver.get_compounds",
        fake_get_compounds,
        raising=True,
    )

    logged_warnings = []

    class FakeLogger:
        @staticmethod
        def warning(msg):
            logged_warnings.append(msg)

    monkeypatch.setattr(
        "cholla_chem.resolvers.pubchem_resolver.pubchem_resolver.logger",
        FakeLogger,
        raising=True,
    )

    compound_names = ["first", "second"]
    result = pcr.name_to_smiles_pubchem(compound_names)

    # Because of zip, only the first name will have an entry
    assert result == {}

    # A warning about mismatching lengths should have been logged
    assert any("Mismatching lengths" in msg for msg in logged_warnings)


def test_parse_cid_file_maps_valid_and_invalid_to_dummy():
    content = "ethanol\t702\nwater\tBAD\nacetone\n"
    out = pcr.parse_cid_file(content)

    assert out["ethanol"] == "702"
    assert out["water"] == pcr.DUMMY_CID
    assert out["acetone"] == pcr.DUMMY_CID


def test_create_batch_cid_request_xml_contains_synonyms_for_name():
    xml = pcr.create_batch_cid_request_xml(["ethanol", "water"], "name")

    assert "<PCT-QueryUids_synonyms>" in xml
    assert "<PCT-QueryUids_synonyms_E>ethanol</PCT-QueryUids_synonyms_E>" in xml
    assert "<PCT-QueryUids_synonyms_E>water</PCT-QueryUids_synonyms_E>" in xml
    assert 'PCT-QueryIDExchange_output-type value="cid"' in xml
    assert 'PCT-QueryIDExchange_output-method value="file-pair"' in xml


def test_poll_request_status_success_returns_download_url(monkeypatch):
    success_xml = """
    <PCT-Data>
      <PCT-Status value="success"/>
      <PCT-Download-URL>
        <PCT-Download-URL_url>ftp://example.com/file.txt</PCT-Download-URL_url>
      </PCT-Download-URL>
    </PCT-Data>
    """.strip()

    monkeypatch.setattr(pcr, "send_xml_to_pug", lambda xml: success_xml, raising=True)

    # Ensure we don't actually sleep in the test
    monkeypatch.setattr(pcr.time, "sleep", lambda _: None, raising=True)

    url = pcr.poll_request_status("REQID", check_interval=1, timeout=5)
    assert url == "ftp://example.com/file.txt"


def test_poll_request_status_error_logs_warning_and_returns_none(monkeypatch):
    error_xml = """
    <PCT-Data>
      <PCT-Status value="error"/>
      <PCT-Status-Message>
        <PCT-Status-Message_message>bad request</PCT-Status-Message_message>
      </PCT-Status-Message>
    </PCT-Data>
    """.strip()

    monkeypatch.setattr(pcr, "send_xml_to_pug", lambda xml: error_xml, raising=True)

    warnings = []

    class FakeLogger:
        @staticmethod
        def warning(msg):
            warnings.append(msg)

    monkeypatch.setattr(pcr, "logger", FakeLogger, raising=True)
    monkeypatch.setattr(pcr.time, "sleep", lambda _: None, raising=True)

    url = pcr.poll_request_status("REQID", check_interval=1, timeout=5)
    assert url is None
    assert any("Error retrieving CIDs from PubChem" in w for w in warnings)


def test_download_file_from_pug_converts_ftp_to_https(monkeypatch):
    captured = {}

    def fake_urlopen(req, timeout=None, context=None):
        # `req` is a urllib.request.Request
        captured["url"] = req.full_url
        return _FakeHTTPResponse(b"hello")

    monkeypatch.setattr(pcr, "urlopen", fake_urlopen, raising=True)

    out = pcr.download_file_from_pug("ftp://pubchem.ncbi.nlm.nih.gov/somefile")
    assert out == "hello"
    assert captured["url"].startswith("https://")


def test_get_json_returns_dict_when_get_returns_json_bytes(monkeypatch):
    def fake_get(
        identifier, namespace, domain, operation, output, searchtype, **kwargs
    ):
        return b'{"PC_Compounds": []}'

    monkeypatch.setattr(pcr, "get", fake_get, raising=True)

    out = pcr.get_json(identifier="123", namespace="cid")
    assert out == {"PC_Compounds": []}


def test_get_json_returns_none_on_exception(monkeypatch):
    def fake_get(*args, **kwargs):
        raise RuntimeError("boom")

    monkeypatch.setattr(pcr, "get", fake_get, raising=True)

    out = pcr.get_json(identifier="123", namespace="cid")
    assert out is None


def test_get_compounds_extracts_smiles_and_replaces_dummy_smiles(monkeypatch):
    # Includes:
    # - one real SMILES
    # - one dummy SMILES that should be converted to ""
    # - one compound with no props -> ""
    payload = {
        "PC_Compounds": [
            {"props": [{"urn": {"label": "SMILES"}, "value": {"sval": "CCO"}}]},
            {
                "props": [
                    {"urn": {"label": "SMILES"}, "value": {"sval": pcr.DUMMY_SMILES}}
                ]
            },
            {"props": []},
        ]
    }

    monkeypatch.setattr(pcr, "get_json", lambda **kwargs: payload, raising=True)

    out = pcr.get_compounds(["a", "b", "c"], namespace="name")
    assert out == ["CCO", "", ""]


def test_get_compounds_raises_if_get_json_none(monkeypatch):
    monkeypatch.setattr(pcr, "get_json", lambda **kwargs: None, raising=True)

    with pytest.raises(ValueError, match="Error retrieving compounds from PubChem"):
        pcr.get_compounds("ethanol", namespace="name")


def test_get_compounds_raises_if_missing_pc_compounds(monkeypatch):
    monkeypatch.setattr(pcr, "get_json", lambda **kwargs: {"foo": "bar"}, raising=True)

    with pytest.raises(ValueError, match="Error retrieving compounds from PubChem"):
        pcr.get_compounds("ethanol", namespace="name")


def test_get_async_listkey_flow_polls_until_done(monkeypatch):
    # This tests get() behavior when initial request returns {"Waiting":{"ListKey":...}}
    call_count = {"n": 0}

    def fake_request(
        identifier,
        namespace="cid",
        domain="compound",
        operation=None,
        output="JSON",
        searchtype=None,
        **kwargs,
    ):
        call_count["n"] += 1

        # First call: returns "Waiting" with a ListKey
        if call_count["n"] == 1:
            return _FakeHTTPResponse(b'{"Waiting":{"ListKey":"LK123"}}')

        # Second call: still waiting
        if call_count["n"] == 2:
            return _FakeHTTPResponse(b'{"Waiting":{"ListKey":"LK123"}}')

        # Third call: completed with final data
        return _FakeHTTPResponse(b'{"PC_Compounds": []}')

    monkeypatch.setattr(pcr, "request", fake_request, raising=True)
    monkeypatch.setattr(pcr.time, "sleep", lambda _: None, raising=True)

    out = pcr.get(
        identifier="ethanol",
        namespace="formula",
        domain="compound",
        output="JSON",
        searchtype=None,
    )
    assert json.loads(out.decode()) == {"PC_Compounds": []}
    assert call_count["n"] == 3


def test_batch_retrieve_cids_happy_path(monkeypatch):
    # 1) send_xml_to_pug returns a "waiting" response with reqid
    initial_xml = """
    <PCT-Data>
      <PCT-Waiting>
        <PCT-Waiting_reqid>REQ123</PCT-Waiting_reqid>
      </PCT-Waiting>
    </PCT-Data>
    """.strip()

    monkeypatch.setattr(
        pcr, "create_batch_cid_request_xml", lambda chunk, ns: "<xml/>", raising=True
    )
    monkeypatch.setattr(pcr, "send_xml_to_pug", lambda xml: initial_xml, raising=True)
    monkeypatch.setattr(
        pcr,
        "poll_request_status",
        lambda req_id, check_interval, timeout: "ftp://example.com/file",
        raising=True,
    )
    monkeypatch.setattr(
        pcr,
        "download_file_from_pug",
        lambda url: "ethanol\t702\nwater\t962\n",
        raising=True,
    )

    out = pcr.batch_retrieve_cids(
        ["ethanol", "water"], namespace="name", chunk_size=1000
    )
    assert out == ["702", "962"]


def test_batch_retrieve_cids_handles_missing_reqid_by_returning_dummy(monkeypatch):
    # No reqid + no error details -> should dummy-fill the chunk
    initial_xml = "<PCT-Data></PCT-Data>"

    monkeypatch.setattr(
        pcr, "create_batch_cid_request_xml", lambda chunk, ns: "<xml/>", raising=True
    )
    monkeypatch.setattr(pcr, "send_xml_to_pug", lambda xml: initial_xml, raising=True)

    out = pcr.batch_retrieve_cids(["a", "b", "c"], namespace="name", chunk_size=1000)
    assert out == [pcr.DUMMY_CID, pcr.DUMMY_CID, pcr.DUMMY_CID]


def test_request_with_list_identifier_in_name_namespace_uses_batch_retrieve_and_switches_to_cid(
    monkeypatch,
):
    # We don't want to actually urlopen; just capture what request() passes to it.
    monkeypatch.setattr(
        pcr, "batch_retrieve_cids", lambda ids, namespace: ["1", "2"], raising=True
    )

    captured = {}

    def fake_urlopen(apiurl, postdata=None, context=None):
        captured["apiurl"] = apiurl
        captured["postdata"] = postdata
        return _FakeHTTPResponse(b"{}")

    monkeypatch.setattr(pcr, "urlopen", fake_urlopen, raising=True)
