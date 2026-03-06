from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Set

# Two-letter elements (must check before single-letter)
TWO_LETTER_ELEMENTS: Set[str] = {
    "He",
    "Li",
    "Be",
    "Ne",
    "Na",
    "Mg",
    "Al",
    "Si",
    "Cl",
    "Ar",
    "Ca",
    "Sc",
    "Ti",
    "Cr",
    "Mn",
    "Fe",
    "Co",
    "Ni",
    "Cu",
    "Zn",
    "Ga",
    "Ge",
    "As",
    "Se",
    "Br",
    "Kr",
    "Rb",
    "Sr",
    "Zr",
    "Nb",
    "Mo",
    "Tc",
    "Ru",
    "Rh",
    "Pd",
    "Ag",
    "Cd",
    "In",
    "Sn",
    "Sb",
    "Te",
    "Xe",
    "Cs",
    "Ba",
    "La",
    "Ce",
    "Pr",
    "Nd",
    "Pm",
    "Sm",
    "Eu",
    "Gd",
    "Tb",
    "Dy",
    "Ho",
    "Er",
    "Tm",
    "Yb",
    "Lu",
    "Hf",
    "Ta",
    "Re",
    "Os",
    "Ir",
    "Pt",
    "Au",
    "Hg",
    "Tl",
    "Pb",
    "Bi",
    "Po",
    "At",
    "Rn",
    "Fr",
    "Ra",
    "Ac",
    "Th",
    "Pa",
    "Np",
    "Pu",
    "Am",
    "Cm",
    "Bk",
    "Cf",
}

# Single-letter elements
ONE_LETTER_ELEMENTS: Set[str] = {
    "H",
    "B",
    "C",
    "N",
    "O",
    "F",
    "P",
    "S",
    "K",
    "V",
    "Y",
    "I",
    "W",
    "U",
}

# Standard valences for validation
STANDARD_VALENCES: Dict[str, List[int]] = {
    "H": [1],
    "B": [3],
    "C": [4],
    "N": [3, 5],
    "O": [2],
    "F": [1],
    "Si": [4],
    "P": [3, 5],
    "S": [2, 4, 6],
    "Cl": [1],
    "Br": [1],
    "I": [1],
}


# =============================================================================
# KNOWN FRAGMENTS - Extensible fragment definitions for common structural units
# =============================================================================


@dataclass
class FragmentDefinition:
    """
    Defines a known molecular fragment that can be recognized in structural formulas.

    Attributes:
        pattern: The text pattern to match (e.g., 'C6H5')
        smiles: The SMILES representation of the fragment with attachment point(s)
        attachment_points: List of atom indices (0-indexed) where bonds can form
        is_aromatic: Whether this fragment contains aromatic atoms
        description: Human-readable description
    """

    pattern: str
    smiles: str
    attachment_points: List[int]  # Indices in the SMILES where attachments occur
    is_aromatic: bool = False
    description: str = ""


KNOWN_FRAGMENTS: Dict[str, FragmentDefinition] = {
    # Phenyl group (benzene with one H removed)
    "C6H5": FragmentDefinition(
        pattern="C6H5",
        smiles="c1ccccc1",
        attachment_points=[0],  # Attach at first carbon
        is_aromatic=True,
        description="phenyl",
    ),
    # Cyclopentadienyl
    "C5H5": FragmentDefinition(
        pattern="C5H5",
        smiles="c1cccc1",
        attachment_points=[0],
        is_aromatic=True,
        description="cyclopentadienyl",
    ),
    # Cyclohexyl
    "C6H11": FragmentDefinition(
        pattern="C6H11",
        smiles="C1CCCCC1",
        attachment_points=[0],
        is_aromatic=False,
        description="cyclohexyl",
    ),
    # Cyclopentyl
    "C5H9": FragmentDefinition(
        pattern="C5H9",
        smiles="C1CCCC1",
        attachment_points=[0],
        is_aromatic=False,
        description="cyclopentyl",
    ),
    # Naphthyl (1-position)
    "C10H7": FragmentDefinition(
        pattern="C10H7",
        smiles="c1ccc2ccccc2c1",
        attachment_points=[0],
        is_aromatic=True,
        description="1-naphthyl",
    ),
    # Carboxylic acid
    "COOH": FragmentDefinition(
        pattern="COOH",
        smiles="C(=O)O",
        attachment_points=[0],
        is_aromatic=False,
        description="carboxylic",
    ),
    # Carboxylic acid
    "HOOC": FragmentDefinition(
        pattern="HOOC",
        smiles="C(=O)O",
        attachment_points=[0],
        is_aromatic=False,
        description="carboxylic",
    ),
    # Primary Amide
    "CONH2": FragmentDefinition(
        pattern="CONH2",
        smiles="",
        attachment_points=[0],
        is_aromatic=False,
        description="primary amide",
    ),
    # Methyl ester
    "COOCH3": FragmentDefinition(
        pattern="COOCH3",
        smiles="",
        attachment_points=[0],
        is_aromatic=False,
        description="methyl ester",
    ),
    # Aldehyde
    "COCH3": FragmentDefinition(
        pattern="COCH3",
        smiles="",
        attachment_points=[0],
        is_aromatic=False,
        description="aldehyde",
    ),
    # Non-terminal aldehyde
    "COCH2": FragmentDefinition(
        pattern="COCH2",
        smiles="",
        attachment_points=[0, 2],
        is_aromatic=False,
        description="non-terminal aldehyde",
    ),
    # Ethyl ester
    "COOCH2CH3": FragmentDefinition(
        pattern="COOCH2CH3",
        smiles="",
        attachment_points=[0],
        is_aromatic=False,
        description="ethyl ester",
    ),
    # Acyl Chloride
    "COCl": FragmentDefinition(
        pattern="COCl",
        smiles="",
        attachment_points=[0],
        is_aromatic=False,
        description="acyl chloride",
    ),
    # Acyl bromide
    "COBr": FragmentDefinition(
        pattern="COBr",
        smiles="",
        attachment_points=[0],
        is_aromatic=False,
        description="acyl bromide",
    ),
    #
    "C2H5": FragmentDefinition(
        pattern="C2H5",
        smiles="",
        attachment_points=[0],
        is_aromatic=False,
        description="",
    ),
    #
    "C3H7": FragmentDefinition(
        pattern="C2H5",
        smiles="",
        attachment_points=[0],
        is_aromatic=False,
        description="",
    ),
    #
    "C4H9": FragmentDefinition(
        pattern="C2H5",
        smiles="",
        attachment_points=[0],
        is_aromatic=False,
        description="",
    ),
    "NO2": FragmentDefinition(
        pattern="NO2",
        smiles="",
        attachment_points=[0],
        is_aromatic=False,
        description="",
    ),
    "O2N": FragmentDefinition(
        pattern="O2N",
        smiles="",
        attachment_points=[2],
        is_aromatic=False,
        description="",
    ),
}
