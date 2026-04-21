"""cholla_chem initialization."""

from importlib.metadata import PackageNotFoundError, version

from cholla_chem.main import (
    ChemicalNameResolver,
    ChemSpiPyResolver,
    CIRpyNameResolver,
    InorganicShorthandNameResolver,
    ManualNameResolver,
    OpsinNameResolver,
    PubChemNameResolver,
    PubChemNameResolverBatch,
    StructuralFormulaNameResolver,
    resolve_compounds_to_smiles,
)
from cholla_chem.name_manipulation.name_correction.dataclasses import CorrectorConfig
from cholla_chem.name_manipulation.name_correction.name_corrector import (
    ChemNameCorrector,
)

__all__ = [
    "resolve_compounds_to_smiles",
    "ChemSpiPyResolver",
    "ChemicalNameResolver",
    "ManualNameResolver",
    "OpsinNameResolver",
    "PubChemNameResolverBatch",
    "PubChemNameResolver",
    "StructuralFormulaNameResolver",
    "CIRpyNameResolver",
    "InorganicShorthandNameResolver",
    "ChemNameCorrector",
    "CorrectorConfig",
]

try:
    __version__ = version("cholla_chem")
except PackageNotFoundError:  # pragma: no cover
    __version__ = "unknown"
