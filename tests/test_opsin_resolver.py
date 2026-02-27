import multiprocessing
import os
import subprocess
import sys

import pytest

PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

from cholla_chem.resolvers.opsin_resolver.opsin_resolver import (  # noqa: E402
    OpsinResult,
    name_to_smiles_opsin,
    run_opsin,
)


def _has_java() -> bool:
    try:
        p = subprocess.run(
            ["java", "-version"], stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        return p.returncode == 0
    except OSError:
        return False


def _f(b):
    return run_opsin(b[0])


CHEMICAL_NAMES = (
    "ethane",
    "water",
    "phenylalanine",
    "(2S)-2-azaniumyl-3-phenylpropanoate",
    "l-mercapto-Z-thiapropane",
)

CHEMICAL_SMILES = (
    "CC",
    "O",
    "N[C@@H](CC1=CC=CC=C1)C(=O)O",
    "[NH3+][C@H](C(=O)[O-])CC1=CC=CC=C1",
    "",
)

CHEMICAL_STDINCHIS = (
    "InChI=1S/C2H6/c1-2/h1-2H3",
    "InChI=1S/H2O/h1H2",
    "InChI=1S/C9H11NO2/c10-8(9(11)12)6-7-4-2-1-3-5-7/h1-5,8H,6,10H2,(H,11,12)/t8-/m0/s1",
    "InChI=1S/C9H11NO2/c10-8(9(11)12)6-7-4-2-1-3-5-7/h1-5,8H,6,10H2,(H,11,12)/t8-/m0/s1",
    "",
)

CHEMICAL_STDINCHIKEYS = (
    "OTMSDBZUPAUEDD-UHFFFAOYSA-N",
    "XLYOFNOQVPJJNP-UHFFFAOYSA-N",
    "COLNVLDHVKWLRT-QMMMGPOBSA-N",
    "COLNVLDHVKWLRT-QMMMGPOBSA-N",
    "",
)

CHEMICAL_INCHI_FIXEDH = (
    "InChI=1/C2H6/c1-2/h1-2H3",
    "InChI=1/H2O/h1H2",
    "InChI=1/C9H11NO2/c10-8(9(11)12)6-7-4-2-1-3-5-7/h1-5,8H,6,10H2,(H,11,12)/t8-/m0/s1/f/h11H",
    "InChI=1/C9H11NO2/c10-8(9(11)12)6-7-4-2-1-3-5-7/h1-5,8H,6,10H2,(H,11,12)/t8-/m0/s1/f/h10H",
    "",
)

CHEMICAL_EXTENDEDSMILES = (
    "CC |$_AV:1;2$|",
    "O |$_AV:O$|",
    "N[C@@H](CC1=CC=CC=C1)C(=O)O |$_AV:N;alpha;beta;1;2;3;4;5;6;;O';O$|",
    "[NH3+][C@H](C(=O)[O-])CC1=CC=CC=C1 |$_AV:1;2;1;O';O;3;1;2;3;4;5;6$|",
    "",
)

CHEMICAL_ERRORS = (
    "",
    "",
    "",
    "",
    "l-mercapto-Z-thiapropane is unparsable due to the following being uninterpretable: l-mercapto-Z-thiapropane The following was not parseable: mercapto-Z-thiapropane",
)

CHEMICAL_INFO = [
    {
        "name": name,
        "smiles": smiles,
        "stdinchi": stdinchi,
        "stdinchikey": stdinchikey,
        "inchi_fixedH": inchi_fixedH,
        "extendedsmiles": extendedsmiles,
        "errors": error,
    }
    for name, smiles, stdinchi, stdinchikey, inchi_fixedH, extendedsmiles, error in zip(
        CHEMICAL_NAMES,
        CHEMICAL_SMILES,
        CHEMICAL_STDINCHIS,
        CHEMICAL_STDINCHIKEYS,
        CHEMICAL_INCHI_FIXEDH,
        CHEMICAL_EXTENDEDSMILES,
        CHEMICAL_ERRORS,
    )
]


@pytest.mark.skipif(not _has_java(), reason="Java is not installed or not on PATH")
def test_multiprocessing():
    """py2opsin should safely work when run with multiprocessing"""
    with multiprocessing.Pool(2) as pool:
        res = pool.map(_f, [("methanol", 0), ("ethanol", 1)])
    assert [item.outputs[0] for item in res] == ["CO", "C(C)O"]
    assert [item.errors[0] for item in res] == ["", ""]
    assert [item.returncode for item in res] == [0, 0]


@pytest.mark.skipif(not _has_java(), reason="Java is not installed or not on PATH")
def test_name_to_smiles():
    """
    Tests converting IUPAC names to SMILES strings
    """
    for test_info in CHEMICAL_INFO:
        opsin_smiles = run_opsin(test_info["name"])
        print("OPSIN outputs:", opsin_smiles.outputs)
        print("OPSIN errors:", opsin_smiles.errors)
        print("OPSIN returncode:", opsin_smiles.returncode)

        assert opsin_smiles.outputs[0] == test_info["smiles"]
        assert opsin_smiles.errors[0] == test_info["errors"]
        assert opsin_smiles.returncode == 0

    test_list_smiles = run_opsin(CHEMICAL_NAMES)
    assert test_list_smiles.outputs == list(CHEMICAL_SMILES)
    assert test_list_smiles.errors == list(CHEMICAL_ERRORS)
    assert test_list_smiles.returncode == 0


@pytest.mark.skipif(not _has_java(), reason="Java is not installed or not on PATH")
def test_name_to_extendedsmiles():
    """
    Tests converting IUPAC names to Extended SMILES
    """
    for test_info in CHEMICAL_INFO:
        opsin_smiles = run_opsin(test_info["name"], output_format="ExtendedSMILES")
        assert opsin_smiles.outputs[0] == test_info["extendedsmiles"]
        assert opsin_smiles.errors[0] == test_info["errors"]
        assert opsin_smiles.returncode == 0

    test_list_extendedsmi = run_opsin(CHEMICAL_NAMES, output_format="ExtendedSMILES")
    assert test_list_extendedsmi.outputs == list(CHEMICAL_EXTENDEDSMILES)
    assert test_list_extendedsmi.errors == list(CHEMICAL_ERRORS)
    assert test_list_extendedsmi.returncode == 0


@pytest.mark.skipif(not _has_java(), reason="Java is not installed or not on PATH")
def test_name_to_stdinchi():
    """
    Tests converting IUPAC names to standard InChI
    """
    for test_info in CHEMICAL_INFO:
        opsin_smiles = run_opsin(test_info["name"], output_format="StdInChI")
        assert opsin_smiles.outputs[0] == test_info["stdinchi"]
        assert opsin_smiles.errors[0] == test_info["errors"]
        assert opsin_smiles.returncode == 0

    test_list_stdinchis = run_opsin(CHEMICAL_NAMES, output_format="StdInChI")
    assert test_list_stdinchis.outputs == list(CHEMICAL_STDINCHIS)
    assert test_list_stdinchis.errors == list(CHEMICAL_ERRORS)
    assert test_list_stdinchis.returncode == 0


@pytest.mark.skipif(not _has_java(), reason="Java is not installed or not on PATH")
def test_name_to_stdinchikey():
    """
    Tests converting IUPAC names to standard InChI keys
    """
    for test_info in CHEMICAL_INFO:
        opsin_smiles = run_opsin(test_info["name"], output_format="StdInChIKey")
        assert opsin_smiles.outputs[0] == test_info["stdinchikey"]
        assert opsin_smiles.errors[0] == test_info["errors"]
        assert opsin_smiles.returncode == 0

    test_list_stdinchikeys = run_opsin(CHEMICAL_NAMES, output_format="StdInChIKey")
    assert test_list_stdinchikeys.outputs == list(CHEMICAL_STDINCHIKEYS)
    assert test_list_stdinchikeys.errors == list(CHEMICAL_ERRORS)
    assert test_list_stdinchikeys.returncode == 0


@pytest.mark.skipif(not _has_java(), reason="Java is not installed or not on PATH")
def test_name_to_inchi_fixedh():
    """
    Tests converting IUPAC names to standard InChI with fixed H
    """
    for test_info in CHEMICAL_INFO:
        opsin_inchi = run_opsin(test_info["name"], output_format="InChI")
        assert opsin_inchi.outputs[0] == test_info["inchi_fixedH"]
        assert opsin_inchi.errors[0] == test_info["errors"]
        assert opsin_inchi.returncode == 0

    test_list_inchi = run_opsin(CHEMICAL_NAMES, output_format="InChI")
    assert test_list_inchi.outputs == list(CHEMICAL_INCHI_FIXEDH)
    assert test_list_inchi.errors == list(CHEMICAL_ERRORS)
    assert test_list_inchi.returncode == 0


@pytest.mark.skipif(not _has_java(), reason="Java is not installed or not on PATH")
def test_allow_multiple_options():
    """
    Test whether run_opsin can handle multiple arguments passed to it
    """
    test_inchi = run_opsin(
        chemical_name="ethane",
        output_format="InChI",
        allow_acid=True,
        allow_radicals=True,
        allow_bad_stereo=True,
        wildcard_radicals=True,
    )

    assert test_inchi.outputs[0] == "InChI=1/C2H6/c1-2/h1-2H3"
    assert test_inchi.errors[0] == ""
    assert test_inchi.returncode == 0


@pytest.mark.skipif(not _has_java(), reason="Java is not installed or not on PATH")
def test_list_with_errors():
    """
    Test whether OPSIN will return a list if there is at least one failed translation
    """
    list_with_errors = ["methane", "ethane", "blah", "water"]
    correct_list = ["C", "CC", "", "O"]
    errors_list = [
        "",
        "",
        "blah is unparsable due to the following being uninterpretable: blah The following was not parseable: blah",
        "",
    ]
    smiles_list = run_opsin(list_with_errors)
    assert smiles_list.outputs == correct_list
    assert smiles_list.errors == errors_list
    assert smiles_list.returncode == 0


def test_name_to_smiles_opsin_success(monkeypatch):
    """name_to_smiles_opsin should map each input name to its SMILES using OPSIN."""

    captured_args = {}

    def fake_run_opsin(
        chemical_name,
        output_format,
        failure_analysis,
        allow_acid,
        allow_radicals,
        allow_bad_stereo,
        wildcard_radicals,
    ):
        captured_args["chemical_name"] = chemical_name
        captured_args["output_format"] = output_format
        captured_args["failure_analysis"] = failure_analysis
        captured_args["allow_acid"] = allow_acid
        captured_args["allow_radicals"] = allow_radicals
        captured_args["allow_bad_stereo"] = allow_bad_stereo
        captured_args["wildcard_radicals"] = wildcard_radicals

        outputs = [f"SMILES_{name}" for name in chemical_name]
        errors = [""] * len(chemical_name)
        return OpsinResult(outputs=outputs, errors=errors, returncode=0)

    monkeypatch.setattr(
        "cholla_chem.resolvers.opsin_resolver.opsin_resolver.run_opsin",
        fake_run_opsin,
        raising=True,
    )

    names = ["ethanol", "water", "acetone"]
    result_smiles, result_failures = name_to_smiles_opsin(
        names,
        allow_acid=True,
        allow_radicals=False,
        allow_bad_stereo=True,
        wildcard_radicals=True,
    )

    assert set(result_smiles.keys()) == set(names)
    for name in names:
        assert result_smiles[name] == f"SMILES_{name}"
    assert result_failures == {}

    assert captured_args["chemical_name"] == names
    assert captured_args["output_format"] == "SMILES"
    assert captured_args["failure_analysis"] is True
    assert captured_args["allow_acid"] is True
    assert captured_args["allow_radicals"] is False
    assert captured_args["allow_bad_stereo"] is True
    assert captured_args["wildcard_radicals"] is True


def test_name_to_smiles_opsin_records_failures(monkeypatch):
    """name_to_smiles_opsin should populate the failure dict when OPSIN returns messages."""

    def fake_run_opsin(chemical_name, **kwargs):
        outputs = [""] * len(chemical_name)
        errors = [f"Error for {name}" for name in chemical_name]
        return OpsinResult(outputs=outputs, errors=errors, returncode=0)

    monkeypatch.setattr(
        "cholla_chem.resolvers.opsin_resolver.opsin_resolver.run_opsin",
        fake_run_opsin,
        raising=True,
    )

    names = ["bad1", "bad2"]
    result_smiles, result_failures = name_to_smiles_opsin(names)

    assert result_smiles == {}
    assert set(result_failures.keys()) == set(names)
    for name in names:
        assert result_failures[name] == f"Error for {name}"


def test_name_to_smiles_opsin_strips_newlines(monkeypatch):
    """name_to_smiles_opsin should strip newline characters before passing to py2opsin."""

    seen_chemical_names = []

    def fake_run_opsin(chemical_name, **kwargs):
        seen_chemical_names.extend(chemical_name)
        outputs = [f"SMILES_{name}" for name in chemical_name]
        errors = [""] * len(chemical_name)
        return OpsinResult(outputs=outputs, errors=errors, returncode=0)

    monkeypatch.setattr(
        "cholla_chem.resolvers.opsin_resolver.opsin_resolver.run_opsin",
        fake_run_opsin,
        raising=True,
    )

    names_with_newlines = ["ethanol\n", "water\n"]
    result_smiles, result_failures = name_to_smiles_opsin(names_with_newlines)

    assert seen_chemical_names == ["ethanol", "water"]
    assert set(result_smiles.keys()) == set(names_with_newlines)
    assert result_failures == {}


def test_name_to_smiles_opsin_mismatched_lengths_logs_warning_and_returns_empty(
    monkeypatch,
):
    """If OPSIN returns mismatched lengths, the function should log a warning and return empty dicts."""

    def fake_run_opsin(chemical_name, **kwargs):
        return OpsinResult(
            outputs=["SMILES_only_one"],
            errors=["err1", "err2"],
            returncode=0,
        )

    monkeypatch.setattr(
        "cholla_chem.resolvers.opsin_resolver.opsin_resolver.run_opsin",
        fake_run_opsin,
        raising=True,
    )

    warnings = []

    class FakeLogger:
        @staticmethod
        def warning(msg):
            warnings.append(msg)

    monkeypatch.setattr(
        "cholla_chem.resolvers.opsin_resolver.opsin_resolver.logger",
        FakeLogger,
        raising=True,
    )

    names = ["n1", "n2"]
    result_smiles, result_failures = name_to_smiles_opsin(names)

    assert result_smiles == {}
    assert result_failures == {}
    assert len(warnings) == 1
    assert "Mismatching lengths" in warnings[0]
