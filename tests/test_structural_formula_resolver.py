import os
import sys

import pytest

# Ensure project root is on sys.path so we can import cholla_chem modules
PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

Chem = pytest.importorskip("rdkit.Chem")

from cholla_chem.resolvers.structural_formula_resolver.structural_formula_resolver import (  # noqa: E402
    StructuralFormulaConverter,
)

TEST_CASES = [
    # (formula, name, expected_smiles)
    ("(CH3)2CHCH2OH", "isobutanol", "CC(C)CO"),
    ("CH3CH2OH", "ethanol", "CCO"),
    ("CH3C(O)CH3", "acetone", "CC(=O)C"),
    ("CH3OH", "methanol", "CO"),
    ("CH4", "methane", "C"),
    ("CH3CH2CH3", "propane", "CCC"),
    ("CH2=CH2", "ethene", "C=C"),
    ("HC#CH", "ethyne", "C#C"),
    ("CH3COOH", "acetic acid", "CC(=O)O"),
    ("(CH3)3COH", "tert-butanol", "CC(C)(C)O"),
    ("CH3NH2", "methylamine", "CN"),
    ("CH3Cl", "chloromethane", "CCl"),
    ("CHCl3", "chloroform", "ClC(Cl)Cl"),
    ("(CH3)2O", "dimethyl ether", "COC"),
    ("CCl4", "carbon tetrachloride", "ClC(Cl)(Cl)Cl"),
    ("(CH3)2NH", "dimethylamine", "CNC"),
    ("(CH3)3N", "trimethylamine", "CN(C)C"),
    ("CH3CHO", "acetaldehyde", "CC=O"),
    ("HCOOH", "formic acid", "C(=O)O"),
    ("CH3CN", "acetonitrile", "CC#N"),
    ("(CH3)3CCH2OH", "neopentyl alcohol", "CC(C)(C)CO"),
    ("(CH3)2CHCH2CH2OH", "isopentyl alcohol", "CC(C)CCO"),
    ("(CH3CH2)3COH", "3-ethyl-3-pentanol", "CCC(O)(CC)CC"),
    # ("(CH3)3COCH3", "methyl tert-butyl ether", "COC(C)(C)C"),
    ("CH3OCH2CH2OCH3", "1,2-dimethoxyethane", "COCCOC"),
    ("(C2H5)2O", "diethyl ether", "CCOCC"),
    ("(CH3CH2)2O", "diethyl ether", "CCOCC"),
    ("HOCH2CH(OH)CH2OH", "glycerol", "OCC(O)CO"),
    ("C6H5CH2OH", "benzyl alcohol", "OCc1ccccc1"),
    ("C6H5NH2", "aniline", "Nc1ccccc1"),
    ("C6H5NO2", "nitrobenzene", "[O-][N+](=O)c1ccccc1"),
    ("C6H5CHO", "benzaldehyde", "O=Cc1ccccc1"),
    ("C6H5COCH3", "acetophenone", "CC(=O)c1ccccc1"),
    ("C6H5OCH3", "anisole", "COc1ccccc1"),
    ("C6H5CH=CH2", "styrene", "C=Cc1ccccc1"),
    ("C6H5C≡CH", "phenylacetylene", "C#Cc1ccccc1"),
    ("C6H5CN", "benzonitrile", "N#Cc1ccccc1"),
    ("C6H5COOH", "benzoic acid", "OC(=O)c1ccccc1"),
    ("C6H5C(CH3)3", "tert-butylbenzene", "CC(C)(C)c1ccccc1"),
    ("(C6H5)2CH2", "diphenylmethane", "c1ccc(Cc2ccccc2)cc1"),
    ("CHCl3", "chloroform", "ClC(Cl)Cl"),
    ("CCl4", "carbon tetrachloride", "ClC(Cl)(Cl)Cl"),
    ("CF3COOH", "trifluoroacetic acid", "OC(=O)C(F)(F)F"),
    ("BrCH2CH2Br", "1,2-dibromoethane", "BrCCBr"),
    ("(CH3)3CBr", "tert-butyl bromide", "CC(C)(C)Br"),
    ("ClCH2COOH", "chloroacetic acid", "OC(=O)CCl"),
    ("CHBr3", "bromoform", "BrC(Br)Br"),
    ("CF2Cl2", "dichlorodifluoromethane", "FC(F)(Cl)Cl"),
    ("CH3CH(OH)COOH", "lactic acid", "CC(O)C(=O)O"),
    ("(CH3)3CCOOH", "pivalic acid", "CC(C)(C)C(=O)O"),
    # ("CH3COOC2H5", "ethyl acetate", "CCOC(C)=O"),
    ("C6H5COOCH3", "methyl benzoate", "COC(=O)c1ccccc1"),
    # ("(CH3CO)2O", "acetic anhydride", "CC(=O)OC(C)=O"),
    ("CH3COCl", "acetyl chloride", "CC(=O)Cl"),
    ("CH3CONH2", "acetamide", "CC(N)=O"),
    ("C6H5CONH2", "benzamide", "NC(=O)c1ccccc1"),
    ("HOOCCH2CH2COOH", "succinic acid", "OC(=O)CCC(=O)O"),
    # ("HOOC(CH2)4COOH", "adipic acid", "OC(=O)CCCCC(=O)O"),
    ("(C2H5)3N", "triethylamine", "CCN(CC)CC"),
    ("(CH3)2NCHO", "N,N-dimethylformamide", "CN(C)C=O"),
    ("HOCH2CH2NH2", "ethanolamine", "NCCO"),
    # ("H2NCH2CH2NH2", "ethylenediamine", "NCCN"),
    ("(HOCH2CH2)3N", "triethanolamine", "OCCN(CCO)CCO"),
    ("CH2=CHCN", "acrylonitrile", "C=CC#N"),
    ("(CH3)2NNH2", "1,1-dimethylhydrazine", "CN(C)N"),
    ("C6H5NHNH2", "phenylhydrazine", "NNc1ccccc1"),
    ("(CH3)2NCH2CH2N(CH3)2", "TMEDA", "CN(C)CCN(C)C"),
    # ("(CH3)2SO", "dimethyl sulfoxide", "CS(C)=O"),
    # ("(CH3)2SO2", "dimethyl sulfone", "CS(C)(=O)=O"),
    ("CH3SH", "methanethiol", "CS"),
    ("CH3SSCH3", "dimethyl disulfide", "CSSC"),
    ("C6H5SH", "thiophenol", "Sc1ccccc1"),
    ("CH3SCH2CH2OH", "2-(methylthio)ethanol", "CSCCO"),
    ("(C2H5)2S", "diethyl sulfide", "CCSCC"),
    ("(CH3)4Si", "tetramethylsilane", "C[Si](C)(C)C"),
    ("(CH3)3SiCl", "trimethylsilyl chloride", "C[Si](C)(C)Cl"),
    ("(CH3)3SiOSi(CH3)3", "hexamethyldisiloxane", "C[Si](C)(C)O[Si](C)(C)C"),
    # ("(C2H5O)4Si", "tetraethyl orthosilicate", "CCO[Si](OCC)(OCC)OCC"),
    ("C6H5Si(CH3)3", "trimethylphenylsilane", "C[Si](C)(C)c1ccccc1"),
    # ("(CH3O)3PO", "trimethyl phosphate", "COP(=O)(OC)OC"),
    # ("(C2H5O)3PO", "triethyl phosphate", "CCOP(=O)(OCC)OCC"),
    ("(C6H5)3P", "triphenylphosphine", "c1ccc(P(c2ccccc2)c2ccccc2)cc1"),
    ("(C6H5)3PO", "triphenylphosphine oxide", "O=P(c1ccccc1)(c1ccccc1)c1ccccc1"),
    ("B(OH)3", "boric acid", "OB(O)O"),
    ("C6H5B(OH)2", "phenylboronic acid", "OB(O)c1ccccc1"),
    # ("(CH3O)3B", "trimethyl borate", "COB(OC)OC"),
    ("(CH3)2CHCHO", "isobutyraldehyde", "CC(C)C=O"),
    ("(CH3)2CHCOCH3", "3-methyl-2-butanone", "CC(C)C(C)=O"),
    # ("CH3COCH2COCH3", "acetylacetone", "CC(=O)CC(C)=O"),
    # ("CH3COCH2CH(CH3)2", "4-methyl-2-pentanone", "CC(C)CC(C)=O"),
    # ("C6H5COCH2CH3", "propiophenone", "CCC(=O)c1ccccc1"),
    # ("C6H5COC6H5", "benzophenone", "O=C(c1ccccc1)c1ccccc1"),
    ("(CH3)2C=CHCOCH3", "mesityl oxide", "CC(C)=CC(C)=O"),
    ("CH2=CHCH2OH", "allyl alcohol", "OCC=C"),
    ("CH3C≡CH", "propyne", "CC#C"),
    ("CH2=C(CH3)COOCH3", "methyl methacrylate", "COC(=O)C(C)=C"),
    ("CH2=CHCOOH", "acrylic acid", "OC(=O)C=C"),
    ("CH2=CHCOOCH3", "methyl acrylate", "COC(=O)C=C"),
    ("(CH3)2C=C(CH3)2", "2,3-dimethyl-2-butene", "CC(C)=C(C)C"),
    ("CH2=CHCH=CH2", "1,3-butadiene", "C=CC=C"),
    # ("H2NCH2COOH", "glycine", "NCC(=O)O"),
    ("CH3CH(NH2)COOH", "alanine", "CC(N)C(=O)O"),
    ("C6H5CH2CH(NH2)COOH", "phenylalanine", "NC(Cc1ccccc1)C(=O)O"),
    # ("HSCH2CH(NH2)COOH", "cysteine", "NC(CS)C(=O)O"),
    ("(CH3)2CHCH(NH2)COOH", "valine", "CC(C)C(N)C(=O)O"),
    ("HOCH2CH2OCH2CH2OH", "diethylene glycol", "OCCOCCO"),
    ("CH3OCH2CH2OCH2CH2OCH3", "diglyme", "COCCOCCOC"),
    ("C6H5CH(OH)CH3", "1-phenylethanol", "CC(O)c1ccccc1"),
    ("(CH3)2C(OH)C≡CH", "2-methyl-3-butyn-2-ol", "CC(C)(O)C#C"),
    ("NCCH2CH2CN", "succinonitrile", "N#CCCC#N"),
    # ("CH3COCH2CH2COCH3", "acetonylacetone", "CC(=O)CCC(C)=O"),
    # ("HOCH2C(CH2OH)3", "pentaerythritol", "OCC(CO)(CO)CO"),
]


def _canonical_smiles(smiles: str) -> str:
    if smiles is None:
        return ""
    smiles = smiles.strip()
    if not smiles:
        return ""
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None, f"Invalid SMILES produced/expected: {smiles!r}"
    return Chem.MolToSmiles(mol, canonical=True)


@pytest.mark.parametrize("formula,name,expected", TEST_CASES, ids=lambda x: str(x))
def test_structural_formula_converter_matches_expected_smiles(formula, name, expected):
    converter = StructuralFormulaConverter(strict_mode=True)
    got = converter.convert(formula)

    assert got != "", (
        "Converter returned empty string\n"
        f"Name: {name}\n"
        f"Formula: {formula}\n"
        f"Errors: {converter.get_errors()}\n"
    )

    assert _canonical_smiles(got) == _canonical_smiles(expected), (
        "Canonical SMILES mismatch\n"
        f"Name: {name}\n"
        f"Formula: {formula}\n"
        f"Got: {got}\n"
        f"Expected: {expected}\n"
        f"Errors: {converter.get_errors()}\n"
    )
