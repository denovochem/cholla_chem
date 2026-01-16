"""name_to_smiles initialization."""

from .main import (
    ChemSpiPyResolver,
    CIRpyNameResolver,
    ChemicalNameResolver,
    ManualNameResolver,
    OpsinNameResolver,
    PubChemNameResolver,
    StructuralFormulaNameResolver,
    resolve_compounds_to_smiles,
)

from .name_manipulation.name_correction.name_corrector import ChemNameCorrector

__all__ = [
    "resolve_compounds_to_smiles",
    "ChemSpiPyResolver",
    "ChemicalNameResolver",
    "ManualNameResolver",
    "OpsinNameResolver",
    "PubChemNameResolver",
    "StructuralFormulaNameResolver",
    "CIRpyNameResolver",
    "ChemNameCorrector",
]

__version__ = "0.0.1"
