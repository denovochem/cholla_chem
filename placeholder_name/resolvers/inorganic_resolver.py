from __future__ import annotations

import re
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from enum import Enum, auto
from typing import Dict, List, Optional, Tuple


# =============================================================================
# DATA STRUCTURES
# =============================================================================


class LigandType(Enum):
    """Classification of ligand charge types."""

    NEUTRAL = auto()
    ANIONIC = auto()
    CATIONIC = auto()


@dataclass
class LigandInfo:
    """
    Complete information about a ligand.

    Attributes:
        smiles: SMILES representation of the ligand
        denticity: Number of coordination sites (atoms that bind to metal)
        charge: Formal charge of the ligand
        aliases: Alternative names/abbreviations for this ligand
        description: Human-readable description
    """

    smiles: str
    denticity: int = 1
    charge: int = 0
    aliases: Tuple[str, ...] = field(default_factory=tuple)
    description: str = ""

    @property
    def ligand_type(self) -> LigandType:
        """Determine ligand type based on formal charge."""
        if self.charge < 0:
            return LigandType.ANIONIC
        elif self.charge > 0:
            return LigandType.CATIONIC
        return LigandType.NEUTRAL


@dataclass
class MetalInfo:
    """
    Information about a transition metal.

    Attributes:
        symbol: Element symbol (e.g., "Ir")
        name: Full element name (e.g., "Iridium")
        common_oxidation_states: List of typical oxidation states
        atomic_number: Atomic number of the element
    """

    symbol: str
    name: str
    common_oxidation_states: Tuple[int, ...]
    atomic_number: int


@dataclass
class ParsedLigand:
    """
    Represents a ligand as parsed from a complex name.

    Attributes:
        name: The ligand identifier as it appears in the name
        count: Number of this ligand coordinated to the metal
        modifiers: Any prefix/suffix modifiers (e.g., "dF(CF3)" in "dF(CF3)ppy")
    """

    name: str
    count: int = 1
    modifiers: List[str] = field(default_factory=list)

    def __repr__(self) -> str:
        return f"ParsedLigand(name='{self.name}', count={self.count})"


@dataclass
class ParsedComplex:
    """
    Complete parsed representation of a coordination complex.

    Attributes:
        metal: Metal symbol (e.g., "Ir", "Rh")
        ligands: List of ParsedLigand objects
        complex_charge: Overall charge of the complex ion
        multiplicity: Number of formula units (e.g., 2 for dimers like [IrCl(cod)]2)
        counter_ions: List of (ion_name, count) tuples
    """

    metal: str
    ligands: List[ParsedLigand]
    complex_charge: int = 0
    multiplicity: int = 1
    counter_ions: List[Tuple[str, int]] = field(default_factory=list)

    def __repr__(self) -> str:
        return (
            f"ParsedComplex(metal='{self.metal}', "
            f"ligands={self.ligands}, "
            f"charge={self.complex_charge}, "
            f"mult={self.multiplicity}, "
            f"counter_ions={self.counter_ions})"
        )


# =============================================================================
# DATABASES / DICTIONARIES
# =============================================================================

# -----------------------------------------------------------------------------
# Ligand Database
# -----------------------------------------------------------------------------

LIGAND_DATABASE: Dict[str, LigandInfo] = {
    # -------------------------------------------------------------------------
    # Monodentate Neutral Ligands
    # -------------------------------------------------------------------------
    "CO": LigandInfo(
        smiles="[C-]#[O+]",
        denticity=1,
        charge=0,
        aliases=("carbonyl",),
        description="Carbonyl ligand",
    ),
    "PPh3": LigandInfo(
        smiles="c1ccc(P(c2ccccc2)c3ccccc3)cc1",
        denticity=1,
        charge=0,
        aliases=("triphenylphosphine", "Ph3P"),
        description="Triphenylphosphine",
    ),
    "py": LigandInfo(
        smiles="c1ccncc1",
        denticity=1,
        charge=0,
        aliases=("pyridine",),
        description="Pyridine",
    ),
    "NH3": LigandInfo(
        smiles="N",
        denticity=1,
        charge=0,
        aliases=("ammonia", "ammine"),
        description="Ammonia/Ammine ligand",
    ),
    "H2O": LigandInfo(
        smiles="O",
        denticity=1,
        charge=0,
        aliases=("water", "aqua", "aquo"),
        description="Water/Aqua ligand",
    ),
    "MeCN": LigandInfo(
        smiles="CC#N",
        denticity=1,
        charge=0,
        aliases=("acetonitrile", "NCMe"),
        description="Acetonitrile",
    ),
    "PMe3": LigandInfo(
        smiles="CP(C)C",
        denticity=1,
        charge=0,
        aliases=("trimethylphosphine",),
        description="Trimethylphosphine",
    ),
    "PEt3": LigandInfo(
        smiles="CCP(CC)CC",
        denticity=1,
        charge=0,
        aliases=("triethylphosphine",),
        description="Triethylphosphine",
    ),
    "PiPr3": LigandInfo(
        smiles="CC(C)P(C(C)C)C(C)C",
        denticity=1,
        charge=0,
        aliases=("triisopropylphosphine", "P(iPr)3"),
        description="Triisopropylphosphine",
    ),
    "PCy3": LigandInfo(
        smiles="C1(CCCCC1)P(C2CCCCC2)C3CCCCC3",
        denticity=1,
        charge=0,
        aliases=("tricyclohexylphosphine",),
        description="Tricyclohexylphosphine",
    ),
    "PtBu3": LigandInfo(
        smiles="CCCCP(CCCC)CCCC",
        denticity=1,
        charge=0,
        aliases=("tri-tert-butylphosphine", "P(tBu)3"),
        description="Tri-tert-butylphosphine",
    ),
    "P(OMe)3": LigandInfo(
        smiles="COP(C)(OC)OC",
        denticity=1,
        charge=0,
        aliases=("trimethylphosphite",),
        description="Trimethyl phosphite",
    ),
    "P(OEt)3": LigandInfo(
        smiles="O(P(OCC)OCC)CC",
        denticity=1,
        charge=0,
        aliases=("triethylphosphite",),
        description="Triethyl phosphite",
    ),
    "P(OPh)3": LigandInfo(
        smiles="O(P(Oc1ccccc1)Oc2ccccc2)c3ccccc3",
        denticity=1,
        charge=0,
        aliases=("triphenylphosphite",),
        description="Triphenyl phosphite",
    ),
    # N-Heterocyclic Carbenes (NHCs) - hugely important in modern catalysis
    "IMes": LigandInfo(
        smiles="Cc1cc(c(c(c1)C)N2C=CN([C]2)c3c(cc(cc3C)C)C)C",
        denticity=1,
        charge=0,
        aliases=("1,3-bis(2,4,6-trimethylphenyl)imidazol-2-ylidene",),
        description="IMes carbene",
    ),
    "IPr": LigandInfo(
        smiles="",
        denticity=1,
        charge=0,
        aliases=("1,3-bis(2,6-diisopropylphenyl)imidazol-2-ylidene",),
        description="IPr carbene",
    ),
    "SIMes": LigandInfo(
        smiles="Cc1cc(C)c(c(C)c1)N2CCN([C]2)c3c(C)cc(C)cc3C",
        denticity=1,
        charge=0,
        aliases=("1,3-bis(2,4,6-trimethylphenyl)imidazolidin-2-ylidene",),
        description="Saturated IMes carbene",
    ),
    "SIPr": LigandInfo(
        smiles="",
        denticity=1,
        charge=0,
        aliases=("1,3-bis(2,6-diisopropylphenyl)imidazolidin-2-ylidene",),
        description="Saturated IPr carbene",
    ),
    "ICy": LigandInfo(
        smiles="",
        denticity=1,
        charge=0,
        aliases=("1,3-dicyclohexylimidazol-2-ylidene",),
        description="ICy carbene",
    ),
    "ItBu": LigandInfo(
        smiles="",
        denticity=1,
        charge=0,
        aliases=("1,3-di-tert-butylimidazol-2-ylidene",),
        description="ItBu carbene",
    ),
    "IMe": LigandInfo(
        smiles="",
        denticity=1,
        charge=0,
        aliases=("1,3-dimethylimidazol-2-ylidene",),
        description="IMe carbene",
    ),
    "IAd": LigandInfo(
        smiles="",
        denticity=1,
        charge=0,
        aliases=("1,3-bis(adamantyl)imidazol-2-ylidene",),
        description="IAd carbene",
    ),
    # Other common neutral donors
    "THF": LigandInfo(
        smiles="C1CCOC1",
        denticity=1,
        charge=0,
        aliases=("tetrahydrofuran",),
        description="Tetrahydrofuran",
    ),
    "Et2O": LigandInfo(
        smiles="CCOCC",
        denticity=1,
        charge=0,
        aliases=("diethylether", "ether"),
        description="Diethyl ether",
    ),
    "DMF": LigandInfo(
        smiles="CN(C)C=O",
        denticity=1,
        charge=0,
        aliases=("dimethylformamide",),
        description="Dimethylformamide",
    ),
    "DMSO": LigandInfo(
        smiles="CS(=O)C",
        denticity=1,
        charge=0,
        aliases=("dimethylsulfoxide",),
        description="Dimethyl sulfoxide",
    ),
    "NMe3": LigandInfo(
        smiles="CN(C)C",
        denticity=1,
        charge=0,
        aliases=("trimethylamine",),
        description="Trimethylamine",
    ),
    "NEt3": LigandInfo(
        smiles="CCN(CC)CC",
        denticity=1,
        charge=0,
        aliases=("triethylamine", "TEA"),
        description="Triethylamine",
    ),
    "DMAP": LigandInfo(
        smiles="Cc1ccccc1N(C)C",
        denticity=1,
        charge=0,
        aliases=("4-dimethylaminopyridine",),
        description="4-Dimethylaminopyridine",
    ),
    # Isocyanides
    "CNtBu": LigandInfo(
        smiles="CC(C)(C)C#N",
        denticity=1,
        charge=0,
        aliases=("tert-butylisocyanide",),
        description="tert-Butyl isocyanide",
    ),
    "CNXyl": LigandInfo(
        smiles="",
        denticity=1,
        charge=0,
        aliases=("2,6-xylylisocyanide", "2,6-dimethylphenylisocyanide"),
        description="2,6-Xylyl isocyanide",
    ),
    "CNPh": LigandInfo(
        smiles="",
        denticity=1,
        charge=0,
        aliases=("phenylisocyanide",),
        description="Phenyl isocyanide",
    ),
    "CNMe": LigandInfo(
        smiles="",
        denticity=1,
        charge=0,
        aliases=("methylisocyanide",),
        description="Methyl isocyanide",
    ),
    # Carbene ligands (non-NHC)
    "CHPh": LigandInfo(
        smiles="",
        denticity=1,
        charge=0,
        aliases=("benzylidene", "phenylcarbene"),
        description="Benzylidene",
    ),
    # Olefins (η²)
    "ethylene": LigandInfo(
        smiles="C=C",
        denticity=1,
        charge=0,
        aliases=("C2H4", "eth"),
        description="Ethylene (η²)",
    ),
    # Dinitrogen and other small molecules
    "N2": LigandInfo(
        smiles="N#N",
        denticity=1,
        charge=0,
        aliases=("dinitrogen",),
        description="Dinitrogen",
    ),
    "NO": LigandInfo(
        smiles="",
        denticity=1,
        charge=0,
        aliases=("nitrosyl",),
        description="Nitrosyl (neutral counting)",
    ),
    "CS": LigandInfo(
        smiles="",
        denticity=1,
        charge=0,
        aliases=("thiocarbonyl",),
        description="Thiocarbonyl",
    ),
    "SO2": LigandInfo(
        smiles="",
        denticity=1,
        charge=0,
        aliases=("sulfurdioxide",),
        description="Sulfur dioxide",
    ),
    "O2": LigandInfo(
        smiles="", denticity=1, charge=0, aliases=("dioxygen",), description="Dioxygen"
    ),
    "H2": LigandInfo(
        smiles="",
        denticity=1,
        charge=0,
        aliases=("dihydrogen",),
        description="Dihydrogen (η²-H2)",
    ),
    # -------------------------------------------------------------------------
    # Monodentate Anionic Ligands
    # -------------------------------------------------------------------------
    "Cl": LigandInfo(
        smiles="[Cl-]",
        denticity=1,
        charge=-1,
        aliases=("chloro", "chloride", "chlorido"),
        description="Chloride",
    ),
    "Br": LigandInfo(
        smiles="[Br-]",
        denticity=1,
        charge=-1,
        aliases=("bromo", "bromide", "bromido"),
        description="Bromide",
    ),
    "I": LigandInfo(
        smiles="[I-]",
        denticity=1,
        charge=-1,
        aliases=("iodo", "iodide", "iodido"),
        description="Iodide",
    ),
    "F": LigandInfo(
        smiles="[F-]",
        denticity=1,
        charge=-1,
        aliases=("fluoro", "fluoride", "fluorido"),
        description="Fluoride",
    ),
    "H": LigandInfo(
        smiles="[H-]",
        denticity=1,
        charge=-1,
        aliases=("hydrido", "hydride"),
        description="Hydride",
    ),
    "CN": LigandInfo(
        smiles="[C-]#N",
        denticity=1,
        charge=-1,
        aliases=("cyano", "cyanide", "cyanido"),
        description="Cyanide",
    ),
    "OAc": LigandInfo(
        smiles="CC([O-])=O",
        denticity=1,
        charge=-1,
        aliases=("acetato", "acetate", "OAc"),
        description="Acetate",
    ),
    "OMe": LigandInfo(
        smiles="[O-]C",
        denticity=1,
        charge=-1,
        aliases=("methoxo", "methoxide"),
        description="Methoxide",
    ),
    # Common anions
    "OH": LigandInfo(
        smiles="[OH-]",
        denticity=1,
        charge=-1,
        aliases=("hydroxo", "hydroxide", "hydroxido"),
        description="Hydroxide",
    ),
    "OEt": LigandInfo(
        smiles="[O-]C(C)C",
        denticity=1,
        charge=-1,
        aliases=("ethoxo", "ethoxide"),
        description="Ethoxide",
    ),
    "OiPr": LigandInfo(
        smiles="[O-]C(C)C",
        denticity=1,
        charge=-1,
        aliases=("isopropoxo", "isopropoxide"),
        description="Isopropoxide",
    ),
    "OtBu": LigandInfo(
        smiles="[O-]C(C)(C)C",
        denticity=1,
        charge=-1,
        aliases=("tert-butoxo", "tert-butoxide"),
        description="tert-Butoxide",
    ),
    "OPh": LigandInfo(
        smiles="[O-]c1ccccc1",
        denticity=1,
        charge=-1,
        aliases=("phenoxo", "phenoxide", "phenolate"),
        description="Phenoxide",
    ),
    # Alkyls
    "Me": LigandInfo(
        smiles="[CH3-]",
        denticity=1,
        charge=-1,
        aliases=("methyl", "CH3"),
        description="Methyl",
    ),
    "Et": LigandInfo(
        smiles="[CH2-]C",
        denticity=1,
        charge=-1,
        aliases=("ethyl", "C2H5"),
        description="Ethyl",
    ),
    "nBu": LigandInfo(
        smiles="[CH2-]CC",
        denticity=1,
        charge=-1,
        aliases=("n-butyl", "butyl"),
        description="n-Butyl",
    ),
    "Ph": LigandInfo(
        smiles="C1=CC=[C-]C=C1",
        denticity=1,
        charge=-1,
        aliases=("phenyl", "C6H5"),
        description="Phenyl",
    ),
    "Bn": LigandInfo(
        smiles="[CH2-]C1=CC=CC=C1", 
        denticity=1, 
        charge=-1, 
        aliases=("benzyl",), 
        description="Benzyl"
    ),
    "vinyl": LigandInfo(
        smiles="", 
        denticity=1, 
        charge=-1, 
        aliases=("ethenyl",), 
        description="Vinyl"
    ),
    "allyl": LigandInfo(
        smiles="",
        denticity=1,
        charge=-1,
        aliases=("η1-allyl", "propenyl"),
        description="Allyl (σ-bound)",
    ),
    "Np": LigandInfo(
        smiles="CC(C)(C)[CH2-]",
        denticity=1,
        charge=-1,
        aliases=("neopentyl",),
        description="Neopentyl",
    ),
    "Mes": LigandInfo(
        smiles="",
        denticity=1,
        charge=-1,
        aliases=("mesityl", "2,4,6-trimethylphenyl"),
        description="Mesityl",
    ),
    # Silyls
    "SiMe3": LigandInfo(
        smiles="",
        denticity=1,
        charge=-1,
        aliases=("trimethylsilyl", "TMS"),
        description="Trimethylsilyl",
    ),
    "SiPh3": LigandInfo(
        smiles="",
        denticity=1,
        charge=-1,
        aliases=("triphenylsilyl",),
        description="Triphenylsilyl",
    ),
    # Amides
    "NMe2": LigandInfo(
        smiles="",
        denticity=1,
        charge=-1,
        aliases=("dimethylamido",),
        description="Dimethylamide",
    ),
    "NEt2": LigandInfo(
        smiles="",
        denticity=1,
        charge=-1,
        aliases=("diethylamido",),
        description="Diethylamide",
    ),
    "NiPr2": LigandInfo(
        smiles="",
        denticity=1,
        charge=-1,
        aliases=("diisopropylamido",),
        description="Diisopropylamide",
    ),
    "NPh2": LigandInfo(
        smiles="",
        denticity=1,
        charge=-1,
        aliases=("diphenylamido",),
        description="Diphenylamide",
    ),
    "NTMS2": LigandInfo(
        smiles="",
        denticity=1,
        charge=-1,
        aliases=("bis(trimethylsilyl)amido", "HMDS", "N(SiMe3)2"),
        description="Bis(trimethylsilyl)amide",
    ),
    "NHPh": LigandInfo(
        smiles="",
        denticity=1,
        charge=-1,
        aliases=("anilido", "phenylamido"),
        description="Anilide",
    ),
    # Other common anionic
    "SCN": LigandInfo(
        smiles="",
        denticity=1,
        charge=-1,
        aliases=("thiocyanato", "thiocyanate"),
        description="Thiocyanate (S-bound)",
    ),
    "NCS": LigandInfo(
        smiles="",
        denticity=1,
        charge=-1,
        aliases=("isothiocyanato", "isothiocyanate"),
        description="Isothiocyanate (N-bound)",
    ),
    "N3": LigandInfo(
        smiles="",
        denticity=1,
        charge=-1,
        aliases=("azido", "azide"),
        description="Azide",
    ),
    "NO2": LigandInfo(
        smiles="",
        denticity=1,
        charge=-1,
        aliases=("nitrito", "nitrite"),
        description="Nitrite (N-bound nitro)",
    ),
    "ONO": LigandInfo(
        smiles="",
        denticity=1,
        charge=-1,
        aliases=("nitrito-O",),
        description="Nitrite (O-bound nitrito)",
    ),
    "SH": LigandInfo(
        smiles="",
        denticity=1,
        charge=-1,
        aliases=("mercapto", "sulfhydryl", "thiolato"),
        description="Hydrosulfide/Thiolate",
    ),
    "SPh": LigandInfo(
        smiles="",
        denticity=1,
        charge=-1,
        aliases=("thiophenolate", "phenylthiolato"),
        description="Thiophenolate",
    ),
    "SMe": LigandInfo(
        smiles="",
        denticity=1,
        charge=-1,
        aliases=("methylthiolate", "methanethiolato"),
        description="Methylthiolate",
    ),
    "StBu": LigandInfo(
        smiles="",
        denticity=1,
        charge=-1,
        aliases=("tert-butylthiolate",),
        description="tert-Butylthiolate",
    ),
    "OCN": LigandInfo(
        smiles="",
        denticity=1,
        charge=-1,
        aliases=("cyanato", "cyanate"),
        description="Cyanate",
    ),
    "NCO": LigandInfo(
        smiles="",
        denticity=1,
        charge=-1,
        aliases=("isocyanato", "isocyanate"),
        description="Isocyanate",
    ),
    # Carboxylates
    "OBz": LigandInfo(
        smiles="",
        denticity=1,
        charge=-1,
        aliases=("benzoato", "benzoate"),
        description="Benzoate",
    ),
    "OPiv": LigandInfo(
        smiles="",
        denticity=1,
        charge=-1,
        aliases=("pivalato", "pivalate", "trimethylacetate"),
        description="Pivalate",
    ),
    "O2CCF3": LigandInfo(
        smiles="",
        denticity=1,
        charge=-1,
        aliases=("trifluoroacetato", "trifluoroacetate", "TFA"),
        description="Trifluoroacetate",
    ),
    "formate": LigandInfo(
        smiles="",
        denticity=1,
        charge=-1,
        aliases=("formato", "HCO2"),
        description="Formate",
    ),
    # Oxo and related
    "O": LigandInfo(
        smiles="",
        denticity=1,
        charge=-2,
        aliases=("oxo", "oxide", "oxido"),
        description="Oxo",
    ),
    "S": LigandInfo(
        smiles="",
        denticity=1,
        charge=-2,
        aliases=("sulfido", "sulfide"),
        description="Sulfido",
    ),
    "Se": LigandInfo(
        smiles="",
        denticity=1,
        charge=-2,
        aliases=("selenido", "selenide"),
        description="Selenido",
    ),
    "Te": LigandInfo(
        smiles="",
        denticity=1,
        charge=-2,
        aliases=("tellurido", "telluride"),
        description="Tellurido",
    ),
    "NR": LigandInfo(
        smiles="",
        denticity=1,
        charge=-2,
        aliases=("imido",),
        description="Imido (generic)",
    ),
    "NAr": LigandInfo(
        smiles="",
        denticity=1,
        charge=-2,
        aliases=("arylimido",),
        description="Aryl imido",
    ),
    "NtBu": LigandInfo(
        smiles="",
        denticity=1,
        charge=-2,
        aliases=("tert-butylimido",),
        description="tert-Butyl imido",
    ),
    "NAd": LigandInfo(
        smiles="",
        denticity=1,
        charge=-2,
        aliases=("adamantylimido",),
        description="Adamantyl imido",
    ),
    # Borohydrides
    "BH4": LigandInfo(
        smiles="",
        denticity=1,
        charge=-1,
        aliases=("borohydride", "tetrahydroborate"),
        description="Borohydride",
    ),
    # -------------------------------------------------------------------------
    # Bidentate Neutral Ligands
    # -------------------------------------------------------------------------
    "bpy": LigandInfo(
        smiles="c1ccc(-c2ccccn2)nc1",
        denticity=2,
        charge=0,
        aliases=("2,2'-bipyridine", "bipyridine", "bipy"),
        description="2,2'-Bipyridine",
    ),
    "dtbbpy": LigandInfo(
        smiles="CC(C)(C)c1ccnc(-c2cc(C(C)(C)C)ccn2)c1",
        denticity=2,
        charge=0,
        aliases=("4,4'-di-tert-butyl-2,2'-bipyridine", "di-tert-butylbipyridine"),
        description="4,4'-Di-tert-butyl-2,2'-bipyridine",
    ),
    "phen": LigandInfo(
        smiles="c1cnc2c(c1)ccc1cccnc12",
        denticity=2,
        charge=0,
        aliases=("1,10-phenanthroline", "phenanthroline"),
        description="1,10-Phenanthroline",
    ),
    "en": LigandInfo(
        smiles="NCCN",
        denticity=2,
        charge=0,
        aliases=("ethylenediamine",),
        description="Ethylenediamine",
    ),
    "cod": LigandInfo(
        smiles="C1=CCCC=CCC1",
        denticity=2,
        charge=0,
        aliases=("1,5-cyclooctadiene", "cyclooctadiene"),
        description="1,5-Cyclooctadiene (η⁴)",
    ),
    "nbd": LigandInfo(
        smiles="C1C2C=CC1C=C2",
        denticity=2,
        charge=0,
        aliases=("norbornadiene", "2,5-norbornadiene"),
        description="Norbornadiene",
    ),
    "dppe": LigandInfo(
        smiles="c1ccc(P(CCP(c2ccccc2)c2ccccc2)c2ccccc2)cc1",
        denticity=2,
        charge=0,
        aliases=("1,2-bis(diphenylphosphino)ethane",),
        description="1,2-Bis(diphenylphosphino)ethane",
    ),
    "dppm": LigandInfo(
        smiles="c1ccc(P(CP(c2ccccc2)c2ccccc2)c2ccccc2)cc1",
        denticity=2,
        charge=0,
        aliases=("bis(diphenylphosphino)methane",),
        description="Bis(diphenylphosphino)methane",
    ),
    # Bipyridines and derivatives
    "4,4'-dmbpy": LigandInfo(
        smiles="CC1=CC(=NC=C1)C2=NC=CC(=C2)C",
        denticity=2,
        charge=0,
        aliases=("4,4'-dimethyl-2,2'-bipyridine", "dmb"),
        description="4,4'-Dimethyl-2,2'-bipyridine",
    ),
    "5,5'-dmbpy": LigandInfo(
        smiles="CC1=CN=C(C=C1)C2=NC=C(C=C2)C",
        denticity=2,
        charge=0,
        aliases=("5,5'-dimethyl-2,2'-bipyridine",),
        description="5,5'-Dimethyl-2,2'-bipyridine",
    ),
    "dCbpy": LigandInfo(
        smiles="C1=CN=C(C=C1C(=O)[O-])C2=NC=CC(=C2)C(=O)[O-]",
        denticity=2,
        charge=0,
        aliases=("4,4'-dicarboxy-2,2'-bipyridine",),
        description="4,4'-Dicarboxy-2,2'-bipyridine",
    ),
    "dCEbpy": LigandInfo(
        smiles="",
        denticity=2,
        charge=0,
        aliases=("4,4'-dicarboxyethyl-2,2'-bipyridine",),
        description="4,4'-Diethoxycarbonyl-2,2'-bipyridine",
    ),
    # Phenanthroline derivatives
    "dmp": LigandInfo(
        smiles="CC(=O)OI1(c2ccccc2C(=O)O1)(OC(=O)C)OC(=O)C",
        denticity=2,
        charge=0,
        aliases=("2,9-dimethyl-1,10-phenanthroline", "neocuproine"),
        description="2,9-Dimethyl-1,10-phenanthroline",
    ),
    "dpp": LigandInfo(
        smiles="C1=CC=C(C=C1)C2=NC3=C(C=CC4=C3N=C(C=C4)C5=CC=CC=C5)C=C2",
        denticity=2,
        charge=0,
        aliases=("2,9-diphenyl-1,10-phenanthroline", "bathocuproine"),
        description="2,9-Diphenyl-1,10-phenanthroline",
    ),
    "tmp": LigandInfo(
        smiles="CC1=CN=C2C(=C1C)C=CC3=C(C(=CN=C32)C)C",
        denticity=2,
        charge=0,
        aliases=("3,4,7,8-tetramethyl-1,10-phenanthroline",),
        description="3,4,7,8-Tetramethyl-1,10-phenanthroline",
    ),
    # Bipyrimidines and related
    "bpm": LigandInfo(
        smiles="C1=CN=C(N=C1)C2=NC=CC=N2",
        denticity=2,
        charge=0,
        aliases=("2,2'-bipyrimidine", "bipyrimidine"),
        description="2,2'-Bipyrimidine",
    ),
    "bpz": LigandInfo(
        smiles="C1=CN=C(C=N1)C2=NC=CN=C2",
        denticity=2,
        charge=0,
        aliases=("2,2'-bipyrazine", "bipyrazine"),
        description="2,2'-Bipyrazine",
    ),
    # Diphosphines (very important in catalysis)
    "dppp": LigandInfo(
        smiles="P(c1ccccc1)(c2ccccc2)CCCP(c3ccccc3)c4ccccc4",
        denticity=2,
        charge=0,
        aliases=("1,3-bis(diphenylphosphino)propane",),
        description="1,3-Bis(diphenylphosphino)propane",
    ),
    "dppb": LigandInfo(
        smiles="C1=CC=C(C=C1)P(CCCCP(C2=CC=CC=C2)C3=CC=CC=C3)C4=CC=CC=C4",
        denticity=2,
        charge=0,
        aliases=("1,4-bis(diphenylphosphino)butane",),
        description="1,4-Bis(diphenylphosphino)butane",
    ),
    "dppf": LigandInfo(
        smiles="c1ccc(cc1)P(c2ccccc2)C34C5[Fe]3678912(C5C6C74)C3C8C9C1(C23)P(c1ccccc1)c1ccccc1",
        denticity=2,
        charge=0,
        aliases=("1,1'-bis(diphenylphosphino)ferrocene",),
        description="1,1'-Bis(diphenylphosphino)ferrocene",
    ),
    "dcpe": LigandInfo(
        smiles="C1CCC(CC1)P(CCP(C2CCCCC2)C3CCCCC3)C4CCCCC4",
        denticity=2,
        charge=0,
        aliases=("1,2-bis(dicyclohexylphosphino)ethane",),
        description="1,2-Bis(dicyclohexylphosphino)ethane",
    ),
    "dcpm": LigandInfo(
        smiles="C1CCC(CC1)P(CP(C2CCCCC2)C3CCCCC3)C4CCCCC4",
        denticity=2,
        charge=0,
        aliases=("bis(dicyclohexylphosphino)methane",),
        description="Bis(dicyclohexylphosphino)methane",
    ),
    "dmpe": LigandInfo(
        smiles="P(C)(C)CCP(C)C",
        denticity=2,
        charge=0,
        aliases=("1,2-bis(dimethylphosphino)ethane",),
        description="1,2-Bis(dimethylphosphino)ethane",
    ),
    "depe": LigandInfo(
        smiles="CCP(CC)CCP(CC)CC",
        denticity=2,
        charge=0,
        aliases=("1,2-bis(diethylphosphino)ethane",),
        description="1,2-Bis(diethylphosphino)ethane",
    ),
    "dippe": LigandInfo(
        smiles="P(C(C)C)(CCP(C(C)C)C(C)C)C(C)C",
        denticity=2,
        charge=0,
        aliases=("1,2-bis(diisopropylphosphino)ethane",),
        description="1,2-Bis(diisopropylphosphino)ethane",
    ),
    "dtbpe": LigandInfo(
        smiles="CC(C)(C)P(CCP(C(C)(C)C)C(C)(C)C)C(C)(C)C",
        denticity=2,
        charge=0,
        aliases=("1,2-bis(di-tert-butylphosphino)ethane",),
        description="1,2-Bis(di-tert-butylphosphino)ethane",
    ),
    # BINAP and related chiral phosphines (critical in asymmetric catalysis)
    "BINAP": LigandInfo(
        smiles="c1ccc(cc1)P(c2ccccc2)c3ccc4ccccc4c3c5c6ccccc6ccc5P(c7ccccc7)c8ccccc8",
        denticity=2,
        charge=0,
        aliases=("2,2'-bis(diphenylphosphino)-1,1'-binaphthyl",),
        description="BINAP",
    ),
    "TolBINAP": LigandInfo(
        smiles="Cc1ccc(cc1)P(c2ccc(C)cc2)c3ccc4ccccc4c3-c5c(ccc6ccccc56)P(c7ccc(C)cc7)c8ccc(C)cc8",
        denticity=2,
        charge=0,
        aliases=("2,2'-bis(di-p-tolylphosphino)-1,1'-binaphthyl",),
        description="Tol-BINAP",
    ),
    "SEGPHOS": LigandInfo(
        smiles="C1OC2=C(O1)C(=C(C=C2)P(C3=CC=CC=C3)C4=CC=CC=C4)C5=C(C=CC6=C5OCO6)P(C7=CC=CC=C7)C8=CC=CC=C8",
        denticity=2,
        charge=0,
        aliases=("5,5'-bis(diphenylphosphino)-4,4'-bi-1,3-benzodioxole",),
        description="SEGPHOS",
    ),
    "DM-SEGPHOS": LigandInfo(
        smiles="Cc1cc(C)cc(c1)P(c2cc(C)cc(C)c2)c3ccc4OCOc4c3-c5c6OCOc6ccc5P(c7cc(C)cc(C)c7)c8cc(C)cc(C)c8",
        denticity=2,
        charge=0,
        aliases=("5,5'-bis(di(3,5-xylyl)phosphino)-4,4'-bi-1,3-benzodioxole",),
        description="DM-SEGPHOS",
    ),
    "DIFLUORPHOS": LigandInfo(
        smiles="", 
        denticity=2, 
        charge=0, 
        aliases=(), 
        description="DIFLUORPHOS"
    ),
    "MeO-BIPHEP": LigandInfo(
        smiles="",
        denticity=2,
        charge=0,
        aliases=("6,6'-dimethoxy-2,2'-bis(diphenylphosphino)-1,1'-biphenyl",),
        description="MeO-BIPHEP",
    ),
    # Josiphos-type ligands
    "Josiphos": LigandInfo(
        smiles="", denticity=2, charge=0, aliases=(), description="Josiphos"
    ),
    # Other chiral ligands
    "CHIRAPHOS": LigandInfo(
        smiles="",
        denticity=2,
        charge=0,
        aliases=("2,3-bis(diphenylphosphino)butane",),
        description="CHIRAPHOS",
    ),
    "DIOP": LigandInfo(
        smiles="",
        denticity=2,
        charge=0,
        aliases=(
            "2,3-O-isopropylidene-2,3-dihydroxy-1,4-bis(diphenylphosphino)butane",
        ),
        description="DIOP",
    ),
    "DuPhos": LigandInfo(
        smiles="",
        denticity=2,
        charge=0,
        aliases=("1,2-bis(phospholano)benzene",),
        description="DuPhos",
    ),
    "BPE": LigandInfo(
        smiles="",
        denticity=2,
        charge=0,
        aliases=("1,2-bis(phospholano)ethane",),
        description="BPE",
    ),
    # Diamines
    "tmeda": LigandInfo(
        smiles="",
        denticity=2,
        charge=0,
        aliases=("N,N,N',N'-tetramethylethylenediamine", "TMEDA"),
        description="TMEDA",
    ),
    "dach": LigandInfo(
        smiles="",
        denticity=2,
        charge=0,
        aliases=("1,2-diaminocyclohexane", "chxn"),
        description="1,2-Diaminocyclohexane",
    ),
    "dpen": LigandInfo(
        smiles="",
        denticity=2,
        charge=0,
        aliases=("1,2-diphenylethylenediamine",),
        description="1,2-Diphenylethylenediamine",
    ),
    "pn": LigandInfo(
        smiles="",
        denticity=2,
        charge=0,
        aliases=("1,2-diaminopropane", "propylenediamine"),
        description="1,2-Diaminopropane",
    ),
    "bn": LigandInfo(
        smiles="",
        denticity=2,
        charge=0,
        aliases=("2,3-diaminobutane", "2,3-butanediamine"),
        description="2,3-Diaminobutane",
    ),
    # Diimine ligands
    "DAB": LigandInfo(
        smiles="",
        denticity=2,
        charge=0,
        aliases=("1,4-diazabutadiene",),
        description="1,4-Diazabutadiene",
    ),
    "Ar-DAB": LigandInfo(
        smiles="",
        denticity=2,
        charge=0,
        aliases=("N,N'-diaryl-1,4-diazabutadiene",),
        description="N,N'-Diaryl-1,4-diazabutadiene",
    ),
    "dpp-BIAN": LigandInfo(
        smiles="",
        denticity=2,
        charge=0,
        aliases=("bis(2,6-diisopropylphenyl)acenaphthenequinonediimine",),
        description="dpp-BIAN",
    ),
    # Mixed P,N donors
    "PHOX": LigandInfo(
        smiles="",
        denticity=2,
        charge=0,
        aliases=("phosphinooxazoline",),
        description="Phosphinooxazoline",
    ),
    # Schiff bases (common in coordination chemistry)
    "salen-H2": LigandInfo(
        smiles="",
        denticity=2,
        charge=0,
        aliases=("N,N'-bis(salicylidene)ethylenediamine",),
        description="Salen (neutral form)",
    ),
    # Other bidentate neutrals
    "dbm": LigandInfo(
        smiles="C1=CC=C(C=C1)C(=O)CC(=O)C2=CC=CC=C2",
        denticity=2,
        charge=0,
        aliases=("dibenzoylmethane",),
        description="Dibenzoylmethane (neutral)",
    ),
    "dme": LigandInfo(
        smiles="COCCOC",
        denticity=2,
        charge=0,
        aliases=("1,2-dimethoxyethane", "glyme"),
        description="1,2-Dimethoxyethane",
    ),
    "diglyme": LigandInfo(
        smiles="COCCOCCOC",
        denticity=2,
        charge=0,
        aliases=("diethylene glycol dimethyl ether",),
        description="Diglyme",
    ),
    "OPPh3": LigandInfo(
        smiles="O=P(c1ccccc1)(c2ccccc2)c3ccccc3",
        denticity=1,
        charge=0,
        aliases=("triphenylphosphine oxide",),
        description="Triphenylphosphine oxide",
    ),
    # -------------------------------------------------------------------------
    # Bidentate Anionic Ligands
    # -------------------------------------------------------------------------
    "acac": LigandInfo(
        smiles="CC(=O)[CH-]C(=O)C",
        denticity=2,
        charge=-1,
        aliases=("acetylacetonate", "acetylacetonato"),
        description="Acetylacetonate",
    ),
    "ppy": LigandInfo(
        smiles="[c-]1ccccc1-c1ccccn1",
        denticity=2,
        charge=-1,
        aliases=("2-phenylpyridine", "phenylpyridine", "phenylpyridinato"),
        description="2-Phenylpyridinate (C^N cyclometalating)",
    ),
    "dfppy": LigandInfo(
        smiles="Fc1cc(F)c([c-]1)-c1ccccn1",
        denticity=2,
        charge=-1,
        aliases=("2-(2,4-difluorophenyl)pyridine",),
        description="2-(2,4-Difluorophenyl)pyridinate",
    ),
    "F2ppy": LigandInfo(
        smiles="Fc1cc(F)c([c-]1)-c1ccccn1",
        denticity=2,
        charge=-1,
        aliases=("difluorophenylpyridine",),
        description="Difluorophenylpyridinate",
    ),
    "pic": LigandInfo(
        smiles="[O-]C(=O)c1ccccn1",
        denticity=2,
        charge=-1,
        aliases=("picolinate", "picolinato"),
        description="Picolinate",
    ),
    # Cyclometalating C^N ligands
    "bzq": LigandInfo(
        smiles="",
        denticity=2,
        charge=-1,
        aliases=("benzo[h]quinoline", "7,8-benzoquinoline"),
        description="Benzo[h]quinolinate",
    ),
    "thpy": LigandInfo(
        smiles="",
        denticity=2,
        charge=-1,
        aliases=("2-(2-thienyl)pyridine", "thienylpyridine"),
        description="2-(2-Thienyl)pyridinate",
    ),
    "piq": LigandInfo(
        smiles="",
        denticity=2,
        charge=-1,
        aliases=("1-phenylisoquinoline", "phenylisoquinoline"),
        description="1-Phenylisoquinolinate",
    ),
    "mppy": LigandInfo(
        smiles="",
        denticity=2,
        charge=-1,
        aliases=("2-(4-methylphenyl)pyridine", "4-methyl-2-phenylpyridine"),
        description="2-(4-Methylphenyl)pyridinate",
    ),
    "btp": LigandInfo(
        smiles="",
        denticity=2,
        charge=-1,
        aliases=("2-benzothienylpyridine",),
        description="2-Benzothienylpyridinate",
    ),
    "pbpy": LigandInfo(
        smiles="",
        denticity=2,
        charge=-1,
        aliases=("6-phenyl-2,2'-bipyridine",),
        description="6-Phenyl-2,2'-bipyridinate",
    ),
    # Beta-diketonates
    "hfac": LigandInfo(
        smiles="",
        denticity=2,
        charge=-1,
        aliases=("hexafluoroacetylacetonate", "1,1,1,5,5,5-hexafluoroacetylacetonate"),
        description="Hexafluoroacetylacetonate",
    ),
    "tfac": LigandInfo(
        smiles="",
        denticity=2,
        charge=-1,
        aliases=("trifluoroacetylacetonate",),
        description="Trifluoroacetylacetonate",
    ),
    "dbm": LigandInfo(
        smiles="",
        denticity=2,
        charge=-1,
        aliases=("dibenzoylmethanate",),
        description="Dibenzoylmethanate",
    ),
    "thd": LigandInfo(
        smiles="",
        denticity=2,
        charge=-1,
        aliases=("2,2,6,6-tetramethyl-3,5-heptanedionate", "tmhd", "dpm"),
        description="2,2,6,6-Tetramethyl-3,5-heptanedionate",
    ),
    "fod": LigandInfo(
        smiles="",
        denticity=2,
        charge=-1,
        aliases=("6,6,7,7,8,8,8-heptafluoro-2,2-dimethyl-3,5-octanedionate",),
        description="FOD",
    ),
    "trop": LigandInfo(
        smiles="",
        denticity=2,
        charge=-1,
        aliases=("tropolonate",),
        description="Tropolonate",
    ),
    # Carboxylates (bridging/chelating)
    "OAc-bi": LigandInfo(
        smiles="",
        denticity=2,
        charge=-1,
        aliases=("acetate-O,O'",),
        description="Acetate (bidentate)",
    ),
    "CO3": LigandInfo(
        smiles="",
        denticity=2,
        charge=-2,
        aliases=("carbonato", "carbonate"),
        description="Carbonate",
    ),
    "SO4": LigandInfo(
        smiles="",
        denticity=2,
        charge=-2,
        aliases=("sulfato", "sulfate"),
        description="Sulfate",
    ),
    # N,N chelates (anionic)
    "pz": LigandInfo(
        smiles="",
        denticity=1,
        charge=-1,
        aliases=("pyrazolato", "pyrazolate"),
        description="Pyrazolate",
    ),
    "pypz": LigandInfo(
        smiles="",
        denticity=2,
        charge=-1,
        aliases=("3-(2-pyridyl)pyrazolate",),
        description="3-(2-Pyridyl)pyrazolate",
    ),
    "indazolato": LigandInfo(
        smiles="",
        denticity=2,
        charge=-1,
        aliases=("indazolate",),
        description="Indazolate",
    ),
    # N,O chelates
    "quinolinolate": LigandInfo(
        smiles="",
        denticity=2,
        charge=-1,
        aliases=("8-hydroxyquinolinate", "oxinate", "Q"),
        description="8-Quinolinolate",
    ),
    "glycinato": LigandInfo(
        smiles="",
        denticity=2,
        charge=-1,
        aliases=("glycinate", "gly"),
        description="Glycinate",
    ),
    "alaninato": LigandInfo(
        smiles="",
        denticity=2,
        charge=-1,
        aliases=("alaninate", "ala"),
        description="Alaninate",
    ),
    "salicylate": LigandInfo(
        smiles="",
        denticity=2,
        charge=-1,
        aliases=("salicylato", "sal"),
        description="Salicylate",
    ),
    "oxalato-mono": LigandInfo(
        smiles="",
        denticity=2,
        charge=-1,
        aliases=(),
        description="Oxalate (monoanionic, monodentate)",
    ),
    # O,O chelates
    "catecholato": LigandInfo(
        smiles="",
        denticity=2,
        charge=-2,
        aliases=("catecholate", "cat"),
        description="Catecholate",
    ),
    "semiquinone": LigandInfo(
        smiles="", denticity=2, charge=-1, aliases=("sq",), description="Semiquinone"
    ),
    # S,S chelates
    "S2CNMe2": LigandInfo(
        smiles="",
        denticity=2,
        charge=-1,
        aliases=("dimethyldithiocarbamate", "Me2dtc"),
        description="Dimethyldithiocarbamate",
    ),
    "S2CNEt2": LigandInfo(
        smiles="",
        denticity=2,
        charge=-1,
        aliases=("diethyldithiocarbamate", "Et2dtc", "dtc"),
        description="Diethyldithiocarbamate",
    ),
    "S2COEt": LigandInfo(
        smiles="",
        denticity=2,
        charge=-1,
        aliases=("ethylxanthate", "xanthate"),
        description="Ethyl xanthate",
    ),
    "S2PPh2": LigandInfo(
        smiles="",
        denticity=2,
        charge=-1,
        aliases=("diphenylphosphinodithioate", "dtp"),
        description="Diphenylphosphinodithioate",
    ),
    "bdt": LigandInfo(
        smiles="",
        denticity=2,
        charge=-2,
        aliases=("1,2-benzenedithiolate", "benzene-1,2-dithiolate"),
        description="1,2-Benzenedithiolate",
    ),
    "mnt": LigandInfo(
        smiles="",
        denticity=2,
        charge=-2,
        aliases=("maleonitriledithiolate", "1,2-dicyanoethylene-1,2-dithiolate"),
        description="Maleonitriledithiolate",
    ),
    "dmit": LigandInfo(
        smiles="",
        denticity=2,
        charge=-2,
        aliases=("2-thioxo-1,3-dithiole-4,5-dithiolate",),
        description="dmit",
    ),
    "edt": LigandInfo(
        smiles="",
        denticity=2,
        charge=-2,
        aliases=("ethane-1,2-dithiolate", "1,2-ethanedithiolate"),
        description="Ethane-1,2-dithiolate",
    ),
    "tdt": LigandInfo(
        smiles="",
        denticity=2,
        charge=-2,
        aliases=("toluene-3,4-dithiolate",),
        description="Toluene-3,4-dithiolate",
    ),
    # -------------------------------------------------------------------------
    # Tridentate Neutral Ligands
    # -------------------------------------------------------------------------
    # Terpyridines
    "tpy": LigandInfo(
        smiles="c1ccnc(c1)c2cccc(n2)c3ccccn3",
        denticity=3,
        charge=0,
        aliases=("2,2':6',2''-terpyridine", "terpyridine", "terpy"),
        description="2,2':6',2''-Terpyridine",
    ),
    "ttpy": LigandInfo(
        smiles="CC1=CC=C(C=C1)C2=CC(=NC(=C2)C3=CC=CC=N3)C4=CC=CC=N4",
        denticity=3,
        charge=0,
        aliases=("4'-p-tolyl-2,2':6',2''-terpyridine",),
        description="4'-p-Tolyl-terpyridine",
    ),
    "tBu3tpy": LigandInfo(
        smiles="CC(C)(C)c1ccnc(c1)-c2cc(cc(n2)-c3cc(ccn3)C(C)(C)C)C(C)(C)C",
        denticity=3,
        charge=0,
        aliases=("4,4',4''-tri-tert-butyl-2,2':6',2''-terpyridine",),
        description="4,4',4''-Tri-tert-butylterpyridine",
    ),
    # Pincer ligands (hugely important class)
    "PNP": LigandInfo(
        smiles="C1=CC=C(C=C1)P(C2=CC=CC=C2)C3=NC(=CC=C3)P(C4=CC=CC=C4)C5=CC=CC=C5",
        denticity=3,
        charge=0,
        aliases=("bis(phosphino)pyridine",),
        description="PNP pincer",
    ),
    "PCP": LigandInfo(
        smiles="",
        denticity=3,
        charge=0,
        aliases=("bis(phosphino)aryl",),
        description="PCP pincer",
    ),
    "NCN": LigandInfo(
        smiles="",
        denticity=3,
        charge=0,
        aliases=("bis(amino)aryl",),
        description="NCN pincer",
    ),
    "SCS": LigandInfo(
        smiles="",
        denticity=3,
        charge=0,
        aliases=("bis(thio)aryl",),
        description="SCS pincer",
    ),
    "CNC": LigandInfo(
        smiles="",
        denticity=3,
        charge=0,
        aliases=("bis(NHC)pyridine",),
        description="CNC pincer (bis-NHC)",
    ),
    # PyBOX
    "PyBOX": LigandInfo(
        smiles="",
        denticity=3,
        charge=0,
        aliases=("pyridinebisoxazoline", "pybox"),
        description="2,6-Bis(oxazolinyl)pyridine",
    ),
    "iPr-PyBOX": LigandInfo(
        smiles="",
        denticity=3,
        charge=0,
        aliases=(),
        description="2,6-Bis(4-isopropyl-2-oxazolinyl)pyridine",
    ),
    "Ph-PyBOX": LigandInfo(
        smiles="",
        denticity=3,
        charge=0,
        aliases=(),
        description="2,6-Bis(4-phenyl-2-oxazolinyl)pyridine",
    ),
    # Triphosphines
    "triphos": LigandInfo(
        smiles="CC(CP(C1=CC=CC=C1)C2=CC=CC=C2)(CP(C3=CC=CC=C3)C4=CC=CC=C4)CP(C5=CC=CC=C5)C6=CC=CC=C6",
        denticity=3,
        charge=0,
        aliases=("1,1,1-tris(diphenylphosphinomethyl)ethane", "MeC(CH2PPh2)3"),
        description="Triphos",
    ),
    # Triamines
    "dien": LigandInfo(
        smiles="NCCNCCN",
        denticity=3,
        charge=0,
        aliases=("diethylenetriamine",),
        description="Diethylenetriamine",
    ),
    "tacn": LigandInfo(
        smiles="C1CNCCNCCN1",
        denticity=3,
        charge=0,
        aliases=("1,4,7-triazacyclononane",),
        description="1,4,7-Triazacyclononane",
    ),
    "Me3tacn": LigandInfo(
        smiles="CN1CCN(CCN(CC1)C)C",
        denticity=3,
        charge=0,
        aliases=("1,4,7-trimethyl-1,4,7-triazacyclononane",),
        description="1,4,7-Trimethyl-1,4,7-triazacyclononane",
    ),
    # Scorpionate-type (neutral)
    "Tp": LigandInfo(
        smiles="",
        denticity=3,
        charge=-1,
        aliases=("hydrotris(pyrazolyl)borate", "trispyrazolylborate"),
        description="Hydrotris(pyrazolyl)borate",
    ),
    "Tp*": LigandInfo(
        smiles="",
        denticity=3,
        charge=-1,
        aliases=("hydrotris(3,5-dimethylpyrazolyl)borate",),
        description="Hydrotris(3,5-dimethylpyrazolyl)borate",
    ),
    "TpiPr2": LigandInfo(
        smiles="",
        denticity=3,
        charge=-1,
        aliases=("hydrotris(3,5-diisopropylpyrazolyl)borate",),
        description="Hydrotris(3,5-diisopropylpyrazolyl)borate",
    ),
    # Other
    "bpa": LigandInfo(
        smiles="",
        denticity=3,
        charge=0,
        aliases=("bis(2-pyridylmethyl)amine",),
        description="Bis(2-pyridylmethyl)amine",
    ),
    "bpea": LigandInfo(
        smiles="",
        denticity=3,
        charge=0,
        aliases=("N,N-bis(2-pyridylmethyl)ethylamine",),
        description="N,N-Bis(2-pyridylmethyl)ethylamine",
    ),
    # -------------------------------------------------------------------------
    # Tridentate Anionic Ligands
    # -------------------------------------------------------------------------
    # Pincer anionic
    "PCP-": LigandInfo(
        smiles="",
        denticity=3,
        charge=-1,
        aliases=(),
        description="PCP pincer (cyclometalated)",
    ),
    "NCN-": LigandInfo(
        smiles="",
        denticity=3,
        charge=-1,
        aliases=(),
        description="NCN pincer (cyclometalated)",
    ),
    "PNP-": LigandInfo(
        smiles="",
        denticity=3,
        charge=-1,
        aliases=(),
        description="PNP pincer (amido form)",
    ),
    # Bis(imino)pyridine
    "PDI": LigandInfo(
        smiles="",
        denticity=3,
        charge=-1,
        aliases=("bis(imino)pyridine", "pyridinediimine"),
        description="Pyridine-2,6-diimine (reduced form)",
    ),
    # Corroles
    "corrole": LigandInfo(
        smiles="",
        denticity=4,
        charge=-3,
        aliases=(),
        description="Corrole (tridentate in some counting)",
    ),
    # -------------------------------------------------------------------------
    # Tetradentate Neutral Ligands
    # -------------------------------------------------------------------------
    # Macrocycles
    "cyclam": LigandInfo(
        smiles="N1CCNCCCNCCNCCC1",
        denticity=4,
        charge=0,
        aliases=("1,4,8,11-tetraazacyclotetradecane",),
        description="Cyclam",
    ),
    "cyclen": LigandInfo(
        smiles="N1CCNCCNCCNCC1",
        denticity=4,
        charge=0,
        aliases=("1,4,7,10-tetraazacyclododecane",),
        description="Cyclen",
    ),
    "Me4cyclam": LigandInfo(
        smiles="CN1CCCN(C)CCN(C)CCCN(C)CC1",
        denticity=4,
        charge=0,
        aliases=("1,4,8,11-tetramethyl-1,4,8,11-tetraazacyclotetradecane",),
        description="Tetramethylcyclam",
    ),
    # Linear tetradentate
    "trien": LigandInfo(
        smiles="NCCNCCNCCN",
        denticity=4,
        charge=0,
        aliases=("triethylenetetramine",),
        description="Triethylenetetramine",
    ),
    # Tetraphosphines
    "PP3": LigandInfo(
        smiles="",
        denticity=4,
        charge=0,
        aliases=("tris(2-(diphenylphosphino)ethyl)phosphine",),
        description="PP3",
    ),
    # -------------------------------------------------------------------------
    # Tetradentate Anionic Ligands
    # -------------------------------------------------------------------------
    # Porphyrins (crucial ligand class)
    "TPP": LigandInfo(
        smiles="",
        denticity=4,
        charge=-2,
        aliases=("tetraphenylporphyrin", "5,10,15,20-tetraphenylporphyrin"),
        description="Tetraphenylporphyrin",
    ),
    "OEP": LigandInfo(
        smiles="",
        denticity=4,
        charge=-2,
        aliases=("octaethylporphyrin", "2,3,7,8,12,13,17,18-octaethylporphyrin"),
        description="Octaethylporphyrin",
    ),
    "TMP": LigandInfo(
        smiles="",
        denticity=4,
        charge=-2,
        aliases=("tetramesitylporphyrin",),
        description="Tetramesitylporphyrin",
    ),
    "por": LigandInfo(
        smiles="",
        denticity=4,
        charge=-2,
        aliases=("porphyrin", "porphyrinato"),
        description="Porphyrin (generic)",
    ),
    "TPFPP": LigandInfo(
        smiles="",
        denticity=4,
        charge=-2,
        aliases=("tetrakis(pentafluorophenyl)porphyrin",),
        description="Tetrakis(pentafluorophenyl)porphyrin",
    ),
    # Phthalocyanines
    "Pc": LigandInfo(
        smiles="",
        denticity=4,
        charge=-2,
        aliases=("phthalocyanine", "phthalocyaninato"),
        description="Phthalocyanine",
    ),
    # Salen-type
    "salen": LigandInfo(
        smiles="",
        denticity=4,
        charge=-2,
        aliases=("N,N'-ethylenebis(salicylideneiminato)",),
        description="Salen",
    ),
    "salphen": LigandInfo(
        smiles="",
        denticity=4,
        charge=-2,
        aliases=("N,N'-phenylenebis(salicylideneiminato)",),
        description="Salphen",
    ),
    "salophen": LigandInfo(
        smiles="", denticity=4, charge=-2, aliases=(), description="Salophen"
    ),
    "salcn": LigandInfo(
        smiles="",
        denticity=4,
        charge=-2,
        aliases=("N,N'-cyclohexanebis(salicylideneiminato)",),
        description="Salen-cyclohexanediamine",
    ),
    "Jacobsen": LigandInfo(
        smiles="",
        denticity=4,
        charge=-2,
        aliases=("Jacobsen's salen",),
        description="Jacobsen's salen ligand",
    ),
    # -------------------------------------------------------------------------
    # Penta/Hexadentate Neutral Ligands
    # -------------------------------------------------------------------------
    # Pentadentate
    "Me6tren": LigandInfo(
        smiles="CN(C)CCN(CCN(C)C)CCN(C)C",
        denticity=4,
        charge=0,
        aliases=("tris(2-dimethylaminoethyl)amine",),
        description="Tris(2-dimethylaminoethyl)amine",
    ),
    "tpa": LigandInfo(
        smiles="c1ccnc(c1)CN(Cc2ccccn2)Cc3ccccn3",
        denticity=4,
        charge=0,
        aliases=("tris(2-pyridylmethyl)amine",),
        description="Tris(2-pyridylmethyl)amine",
    ),
    # Hexadentate
    "EDTA": LigandInfo(
        smiles="OC(=O)CN(CCN(CC(O)=O)CC(O)=O)CC(O)=O",
        denticity=6,
        charge=-4,
        aliases=("ethylenediaminetetraacetate", "edta"),
        description="Ethylenediaminetetraacetate",
    ),
    "DTPA": LigandInfo(
        smiles="C(CN(CC(=O)O)CC(=O)O)N(CCN(CC(=O)O)CC(=O)O)CC(=O)O",
        denticity=8,
        charge=-5,
        aliases=("diethylenetriaminepentaacetate",),
        description="Diethylenetriaminepentaacetate",
    ),
    "tpen": LigandInfo(
        smiles="C1=CC=NC(=C1)CN(CCN(CC2=CC=CC=N2)CC3=CC=CC=N3)CC4=CC=CC=N4",
        denticity=6,
        charge=0,
        aliases=("N,N,N',N'-tetrakis(2-pyridylmethyl)ethylenediamine",),
        description="TPEN",
    ),
    # -------------------------------------------------------------------------
    # η-Bonded Ligands
    # -------------------------------------------------------------------------
    "Cp": LigandInfo(
        smiles="[cH-]1cccc1",
        denticity=5,
        charge=-1,
        aliases=("cyclopentadienyl", "C5H5"),
        description="Cyclopentadienyl (η⁵)",
    ),
    "Cp*": LigandInfo(
        smiles="C[c-]1c(C)c(C)c(C)c1C",
        denticity=5,
        charge=-1,
        aliases=("pentamethylcyclopentadienyl", "C5Me5", "Cpstar"),
        description="Pentamethylcyclopentadienyl (η⁵)",
    ),
    # Arenes (η⁶)
    "benzene": LigandInfo(
        smiles="c1ccccc1",
        denticity=6,
        charge=0,
        aliases=("η6-benzene", "C6H6"),
        description="Benzene (η⁶)",
    ),
    "p-cymene": LigandInfo(
        smiles="c1cc(ccc1C(C)C)C",
        denticity=6,
        charge=0,
        aliases=("η6-p-cymene", "4-isopropyltoluene"),
        description="p-Cymene (η⁶)",
    ),
    "mesitylene": LigandInfo(
        smiles="Cc1cc(cc(c1)C)C",
        denticity=6,
        charge=0,
        aliases=("η6-mesitylene", "1,3,5-trimethylbenzene"),
        description="Mesitylene (η⁶)",
    ),
    "hexamethylbenzene": LigandInfo(
        smiles="c1(c(c(c(c(c1C)C)C)C)C)C",
        denticity=6,
        charge=0,
        aliases=("η6-C6Me6", "HMB"),
        description="Hexamethylbenzene (η⁶)",
    ),
    "toluene": LigandInfo(
        smiles="Cc1ccccc1",
        denticity=6,
        charge=0,
        aliases=("η6-toluene",),
        description="Toluene (η⁶)",
    ),
    "C6H3Me3": LigandInfo(
        smiles="Cc1cc(cc(c1)C)C",
        denticity=6,
        charge=0,
        aliases=("η6-trimethylbenzene",),
        description="Trimethylbenzene (η⁶)",
    ),
    # Allyl (η³)
    "allyl": LigandInfo(
        smiles="",
        denticity=3,
        charge=-1,
        aliases=("η3-allyl", "η3-C3H5"),
        description="Allyl (η³)",
    ),
    "methallyl": LigandInfo(
        smiles="",
        denticity=3,
        charge=-1,
        aliases=("η3-2-methylallyl", "η3-methallyl"),
        description="Methallyl (η³)",
    ),
    "crotyl": LigandInfo(
        smiles="",
        denticity=3,
        charge=-1,
        aliases=("η3-crotyl", "η3-1-methylallyl"),
        description="Crotyl (η³)",
    ),
    "cinnamyl": LigandInfo(
        smiles="",
        denticity=3,
        charge=-1,
        aliases=("η3-cinnamyl", "η3-3-phenylallyl"),
        description="Cinnamyl (η³)",
    ),
    # Tropylium and related
    "C7H7": LigandInfo(
        smiles="",
        denticity=7,
        charge=1,
        aliases=("tropylium", "cycloheptatrienyl"),
        description="Tropylium (η⁷)",
    ),
    "cht": LigandInfo(
        smiles="C1=C\C/C=C\C=C1",
        denticity=6,
        charge=0,
        aliases=("cycloheptatriene", "η6-C7H8"),
        description="Cycloheptatriene (η⁶)",
    ),
    # Indenyl
    "Ind": LigandInfo(
        smiles="",
        denticity=5,
        charge=-1,
        aliases=("indenyl", "η5-indenyl"),
        description="Indenyl (η⁵)",
    ),
    # Fluorenyl
    "Flu": LigandInfo(
        smiles="",
        denticity=5,
        charge=-1,
        aliases=("fluorenyl", "η5-fluorenyl"),
        description="Fluorenyl (η⁵)",
    ),
    # Cyclooctatetraene
    "cot": LigandInfo(
        smiles="",
        denticity=8,
        charge=-2,
        aliases=("cyclooctatetraene", "η8-C8H8", "COT"),
        description="Cyclooctatetraene (η⁸)",
    ),
    # Butadiene
    "bd": LigandInfo(
        smiles="C=CC=C",
        denticity=4,
        charge=0,
        aliases=("butadiene", "η4-butadiene"),
        description="1,3-Butadiene (η⁴)",
    ),
    "isoprene": LigandInfo(
        smiles="CC(=C)C=C",
        denticity=4,
        charge=0,
        aliases=("η4-isoprene", "2-methylbutadiene"),
        description="Isoprene (η⁴)",
    ),
}


# -----------------------------------------------------------------------------
# Counter Ion Database
# -----------------------------------------------------------------------------

COUNTER_ION_DATABASE: Dict[str, LigandInfo] = {
    "PF6": LigandInfo(
        smiles="F[P-](F)(F)(F)(F)F",
        charge=-1,
        aliases=("hexafluorophosphate",),
        description="Hexafluorophosphate",
    ),
    "BF4": LigandInfo(
        smiles="F[B-](F)(F)F",
        charge=-1,
        aliases=("tetrafluoroborate",),
        description="Tetrafluoroborate",
    ),
    "OTf": LigandInfo(
        smiles="[O-]S(=O)(=O)C(F)(F)F",
        charge=-1,
        aliases=("triflate", "trifluoromethanesulfonate", "CF3SO3"),
        description="Triflate",
    ),
    "ClO4": LigandInfo(
        smiles="[O-]Cl(=O)(=O)=O",
        charge=-1,
        aliases=("perchlorate",),
        description="Perchlorate",
    ),
    "SbF6": LigandInfo(
        smiles="F[Sb-](F)(F)(F)(F)F",
        charge=-1,
        aliases=("hexafluoroantimonate",),
        description="Hexafluoroantimonate",
    ),
    "BArF": LigandInfo(
        smiles="FC(F)(F)c1cc([B-](c2cc(C(F)(F)F)cc(C(F)(F)F)c2)"
        "(c2cc(C(F)(F)F)cc(C(F)(F)F)c2)c2cc(C(F)(F)F)cc(C(F)(F)F)c2)"
        "cc(C(F)(F)F)c1",
        charge=-1,
        aliases=("BArF24", "tetrakis(3,5-bis(trifluoromethyl)phenyl)borate"),
        description="Tetrakis(3,5-bis(trifluoromethyl)phenyl)borate",
    ),
    "BAr4": LigandInfo(
        smiles="c1ccc([B-](c2ccccc2)(c2ccccc2)c2ccccc2)cc1",
        charge=-1,
        aliases=("tetraphenylborate", "BPh4"),
        description="Tetraphenylborate",
    ),
    "NO3": LigandInfo(
        smiles="[O-][N+](=O)[O-]",
        charge=-1,
        aliases=("nitrate",),
        description="Nitrate",
    ),
    "Cl": LigandInfo(
        smiles="[Cl-]", charge=-1, aliases=("chloride",), description="Chloride"
    ),
    "Br": LigandInfo(
        smiles="[Br-]", charge=-1, aliases=("bromide",), description="Bromide"
    ),
    "I": LigandInfo(
        smiles="[I-]", charge=-1, aliases=("iodide",), description="Iodide"
    ),
    "BPh4": LigandInfo(
        smiles="[B-](c1ccccc1)(c2ccccc2)(c3ccccc3)c4ccccc4",
        charge=-1,
        aliases=("tetraphenylborate",),
        description="Tetraphenylborate",
    ),
    "Al(OC(CF3)3)4": LigandInfo(
        smiles="",
        charge=-1,
        aliases=("perfluoro-tert-butoxide aluminate",),
        description="Perfluoro-tert-butoxide aluminate",
    ),
    "BAr4F": LigandInfo(
        smiles="",
        charge=-1,
        aliases=("tetrakis(3,5-bis(trifluoromethyl)phenyl)borate", "BArF"),
        description="BArF",
    ),
    "NTf2": LigandInfo(
        smiles="C(F)(F)(F)S(=O)(=O)[N-]S(=O)(=O)C(F)(F)F",
        charge=-1,
        aliases=("bis(trifluoromethylsulfonyl)imide", "TFSI", "bistriflimide"),
        description="Bis(trifluoromethylsulfonyl)imide",
    ),
    "AsF6": LigandInfo(
        smiles="F[As-](F)(F)(F)(F)F",
        charge=-1,
        aliases=("hexafluoroarsenate",),
        description="Hexafluoroarsenate",
    ),
    "B(C6F5)4": LigandInfo(
        smiles="Fc1c(c(F)c(F)c(F)c1F)[B-](c2c(F)c(F)c(F)c(F)c2F)(c3c(F)c(F)c(F)c(F)c3F)c4c(F)c(F)c(F)c(F)c4F",
        charge=-1,
        aliases=("tetrakis(pentafluorophenyl)borate",),
        description="Tetrakis(pentafluorophenyl)borate",
    ),
    "HSO4": LigandInfo(
        smiles="OS(=O)(=O)[O-]",
        charge=-1,
        aliases=("hydrogensulfate", "bisulfate"),
        description="Hydrogen sulfate",
    ),
    "CF3CO2": LigandInfo(
        smiles="FC(F)(F)C(=O)[O-]",
        charge=-1,
        aliases=("trifluoroacetate", "TFA"),
        description="Trifluoroacetate",
    ),
    "MeSO3": LigandInfo(
        smiles="O=S(=O)([O-])C",
        charge=-1,
        aliases=("mesylate", "methanesulfonate"),
        description="Mesylate",
    ),
    "TsO": LigandInfo(
        smiles="CC1=CC=C(C=C1)S(=O)(=O)[O-]",
        charge=-1,
        aliases=("tosylate", "4-toluenesulfonate", "OTs"),
        description="Tosylate",
    ),
    "ReO4": LigandInfo(
        smiles="O=[Re](O)(O)(O)(O)", charge=-1, aliases=("perrhenate",), description="Perrhenate"
    ),
    "IO4": LigandInfo(
        smiles="O[I](O)(O)(O)O", charge=-1, aliases=("periodate",), description="Periodate"
    ),
    "BH4": LigandInfo(
        smiles="[BH4-]",
        charge=-1,
        aliases=("borohydride", "tetrahydroborate"),
        description="Borohydride (as counter ion)",
    ),
    "AlH4": LigandInfo(
        smiles="[AlH4-]",
        charge=-1,
        aliases=("aluminate", "tetrahydroaluminate"),
        description="Tetrahydroaluminate",
    ),
    "AlCl4": LigandInfo(
        smiles="Cl[Al-](Cl)(Cl)Cl",
        charge=-1,
        aliases=("tetrachloroaluminate",),
        description="Tetrachloroaluminate",
    ),
    "FeCl4": LigandInfo(
        smiles="Cl[Fe-](Cl)(Cl)Cl",
        charge=-1,
        aliases=("tetrachloroferrate",),
        description="Tetrachloroferrate",
    ),
    "CuCl2": LigandInfo(
        smiles="",
        charge=-1,
        aliases=("dichlorocuprate",),
        description="Dichlorocuprate",
    ),
    "ZnCl3": LigandInfo(
        smiles="Cl[Zn-](Cl)Cl",
        charge=-1,
        aliases=("trichlorozincate",),
        description="Trichlorozincate",
    ),
    "GaCl4": LigandInfo(
        smiles="",
        charge=-1,
        aliases=("tetrachlorogallate",),
        description="Tetrachlorogallate",
    ),
    # Cationic counterions (for anionic complexes)
    "Na": LigandInfo(
        smiles="[Na+]", charge=1, aliases=("sodium",), description="Sodium"
    ),
    "K": LigandInfo(
        smiles="[K+]", charge=1, aliases=("potassium",), description="Potassium"
    ),
    "Li": LigandInfo(
        smiles="[Li+]", charge=1, aliases=("lithium",), description="Lithium"
    ),
    "Cs": LigandInfo(
        smiles="[Cs+]", charge=1, aliases=("cesium", "caesium"), description="Cesium"
    ),
    "NBu4": LigandInfo(
        smiles="CCCC[N+](CCCC)(CCCC)CCCC",
        charge=1,
        aliases=("tetrabutylammonium", "TBA", "nBu4N"),
        description="Tetrabutylammonium",
    ),
    "NEt4": LigandInfo(
        smiles="CC[N+](CC)(CC)CC",
        charge=1,
        aliases=("tetraethylammonium", "TEA", "Et4N"),
        description="Tetraethylammonium",
    ),
    "NMe4": LigandInfo(
        smiles="C[N+](C)(C)C",
        charge=1,
        aliases=("tetramethylammonium", "TMA", "Me4N"),
        description="Tetramethylammonium",
    ),
    "PPh4": LigandInfo(
        smiles="c1c(cccc1)[P+](c2ccccc2)(c3ccccc3)c4ccccc4",
        charge=1,
        aliases=("tetraphenylphosphonium", "Ph4P"),
        description="Tetraphenylphosphonium",
    ),
    "PPN": LigandInfo(
        smiles="",
        charge=1,
        aliases=("bis(triphenylphosphine)iminium", "Ph3P=N=PPh3"),
        description="Bis(triphenylphosphine)iminium",
    ),
    "Cp2Fe": LigandInfo(
        smiles="", charge=1, aliases=("ferrocenium", "Fc+"), description="Ferrocenium"
    ),
    "Cp2Co": LigandInfo(
        smiles="", charge=1, aliases=("cobaltocenium",), description="Cobaltocenium"
    ),
    "H3O": LigandInfo(
        smiles="", charge=1, aliases=("hydronium", "oxonium"), description="Hydronium"
    ),
    "NH4": LigandInfo(
        smiles="[NH4+]", charge=1, aliases=("ammonium",), description="Ammonium"
    ),
    "pyH": LigandInfo(
        smiles="", charge=1, aliases=("pyridinium",), description="Pyridinium"
    ),
    "DMAH": LigandInfo(
        smiles="",
        charge=1,
        aliases=("dimethylanilinium",),
        description="Dimethylanilinium",
    ),
}


# -----------------------------------------------------------------------------
# Metal Database
# -----------------------------------------------------------------------------

METAL_DATABASE: Dict[str, MetalInfo] = {
    # Group 6
    "Cr": MetalInfo("Cr", "Chromium", (0, 2, 3, 6), 24),
    "Mo": MetalInfo("Mo", "Molybdenum", (0, 2, 4, 6), 42),
    "W": MetalInfo("W", "Tungsten", (0, 2, 4, 6), 74),
    # Group 7
    "Mn": MetalInfo("Mn", "Manganese", (0, 2, 3, 4, 7), 25),
    "Re": MetalInfo("Re", "Rhenium", (0, 1, 3, 5, 7), 75),
    # Group 8
    "Fe": MetalInfo("Fe", "Iron", (0, 2, 3), 26),
    "Ru": MetalInfo("Ru", "Ruthenium", (0, 2, 3, 4), 44),
    "Os": MetalInfo("Os", "Osmium", (0, 2, 3, 4, 6, 8), 76),
    # Group 9
    "Co": MetalInfo("Co", "Cobalt", (0, 2, 3), 27),
    "Rh": MetalInfo("Rh", "Rhodium", (0, 1, 2, 3), 45),
    "Ir": MetalInfo("Ir", "Iridium", (0, 1, 3, 4), 77),
    # Group 10
    "Ni": MetalInfo("Ni", "Nickel", (0, 2), 28),
    "Pd": MetalInfo("Pd", "Palladium", (0, 2, 4), 46),
    "Pt": MetalInfo("Pt", "Platinum", (0, 2, 4), 78),
    # Group 11
    "Cu": MetalInfo("Cu", "Copper", (1, 2), 29),
    "Ag": MetalInfo("Ag", "Silver", (1,), 47),
    "Au": MetalInfo("Au", "Gold", (1, 3), 79),
    # Group 12
    "Zn": MetalInfo("Zn", "Zinc", (2,), 30),
}


# =============================================================================
# PARSER
# =============================================================================


class ParserError(Exception):
    """Custom exception for parsing errors."""

    pass


class ComplexNameParser:
    """
    Parser for inorganic/organometallic complex names.

    Handles names in common formats such as:
        - [Metal(Ligand)n(Ligand2)m]charge
        - [Metal(Ligand)n]multiplicity
        - [Metal(Ligand)n]CounterIon

    Examples:
        >>> parser = ComplexNameParser()
        >>> result = parser.parse("[IrCl(cod)]2")
        >>> print(result.metal)
        'Ir'
    """

    def __init__(
        self,
        ligand_db: Optional[Dict[str, LigandInfo]] = None,
        metal_db: Optional[Dict[str, MetalInfo]] = None,
        counter_ion_db: Optional[Dict[str, LigandInfo]] = None,
    ) -> None:
        """
        Initialize the parser with chemical databases.

        Args:
            ligand_db: Dictionary mapping ligand abbreviations to LigandInfo.
                      Uses LIGAND_DATABASE if None.
            metal_db: Dictionary mapping metal symbols to MetalInfo.
                     Uses METAL_DATABASE if None.
            counter_ion_db: Dictionary mapping counter ion names to LigandInfo.
                           Uses COUNTER_ION_DATABASE if None.
        """
        self.ligand_db = ligand_db if ligand_db is not None else LIGAND_DATABASE
        self.metal_db = metal_db if metal_db is not None else METAL_DATABASE
        self.counter_ion_db = (
            counter_ion_db if counter_ion_db is not None else COUNTER_ION_DATABASE
        )

    def parse(self, name: str) -> ParsedComplex:
        """
        Parse an inorganic complex name into its components.

        Args:
            name: The complex name string (e.g., "[IrCl(cod)]2")

        Returns:
            ParsedComplex object containing all parsed components

        Raises:
            ParserError: If the name cannot be parsed
        """
        original_name = name.strip()
        working_name = original_name

        # Step 1: Extract trailing multiplicity (e.g., "]2" at end)
        multiplicity, working_name = self._extract_multiplicity(working_name)

        # Step 2: Extract counter ions (after the complex brackets)
        counter_ions, working_name = self._extract_counter_ions(working_name)

        # Step 3: Extract complex charge (e.g., "]+" or "]2-")
        complex_charge, working_name = self._extract_charge(working_name)

        # Step 4: Remove outer brackets
        working_name = self._strip_brackets(working_name)

        # Step 5: Extract metal symbol
        metal, ligand_string = self._extract_metal(working_name)

        # Step 6: Parse ligands
        ligands = self._parse_ligand_string(ligand_string)

        return ParsedComplex(
            metal=metal,
            ligands=ligands,
            complex_charge=complex_charge,
            multiplicity=multiplicity,
            counter_ions=counter_ions,
        )

    def _extract_multiplicity(self, name: str) -> Tuple[int, str]:
        """
        Extract multiplicity from end of name (e.g., "]2").

        Args:
            name: Complex name string

        Returns:
            Tuple of (multiplicity, remaining_string)
        """
        # Pattern: ends with ]<number> where number is the multiplicity
        match = re.search(r"\](\d+)$", name)
        if match:
            multiplicity = int(match.group(1))
            # Keep the closing bracket, remove the number
            return multiplicity, name[: match.start() + 1]
        return 1, name

    def _extract_counter_ions(self, name: str) -> Tuple[List[Tuple[str, int]], str]:
        """
        Extract counter ions from after the complex brackets.

        Args:
            name: Complex name string

        Returns:
            Tuple of (list of (ion_name, count) tuples, remaining_string)
        """
        counter_ions: List[Tuple[str, int]] = []

        # Find where the complex ends (after closing bracket)
        if not name.endswith("]"):
            # Look for content after the last ]
            match = re.search(r"\]([A-Za-z0-9()]+)$", name)
            if match:
                counter_string = match.group(1)
                remaining = name[: match.start() + 1]

                # Try to match known counter ions
                for ion_name in sorted(
                    self.counter_ion_db.keys(), key=len, reverse=True
                ):
                    # Look for the ion with optional count
                    ion_pattern = rf"({re.escape(ion_name)})(\d*)"
                    ion_match = re.search(ion_pattern, counter_string)
                    if ion_match:
                        count = int(ion_match.group(2)) if ion_match.group(2) else 1
                        counter_ions.append((ion_name, count))
                        counter_string = counter_string.replace(
                            ion_match.group(0), "", 1
                        )

                return counter_ions, remaining

        return counter_ions, name

    def _extract_charge(self, name: str) -> Tuple[int, str]:
        """
        Extract complex charge notation (e.g., "]+" or "]2-").

        Args:
            name: Complex name string

        Returns:
            Tuple of (charge, remaining_string)
        """
        # Pattern: ]<optional_number><+/->
        match = re.search(r"\](\d*)([+-])$", name)
        if match:
            charge_magnitude = int(match.group(1)) if match.group(1) else 1
            charge_sign = 1 if match.group(2) == "+" else -1
            charge = charge_magnitude * charge_sign
            return charge, name[: match.start() + 1]
        return 0, name

    def _strip_brackets(self, name: str) -> str:
        """
        Remove outer square brackets from complex name.

        Args:
            name: Complex name string

        Returns:
            String with outer brackets removed
        """
        name = name.strip()
        if name.startswith("[") and name.endswith("]"):
            return name[1:-1]
        return name

    def _extract_metal(self, name: str) -> Tuple[str, str]:
        """
        Extract metal symbol from beginning of name.

        Args:
            name: Complex name string (without outer brackets)

        Returns:
            Tuple of (metal_symbol, remaining_ligand_string)

        Raises:
            ParserError: If no known metal is found
        """
        # Try metals sorted by length (longest first) to handle cases like
        # "Ir" vs "I" (iodine)
        for metal_symbol in sorted(self.metal_db.keys(), key=len, reverse=True):
            if name.startswith(metal_symbol):
                return metal_symbol, name[len(metal_symbol) :]

        raise ParserError(f"Could not identify metal in: {name}")

    def _parse_ligand_string(self, ligand_str: str) -> List[ParsedLigand]:
        """
        Parse the ligand portion of a complex name.

        Handles:
            - Parenthesized ligands: (cod), (ppy)2
            - Direct ligands: Cl, Cl2
            - Complex nested names: dF(CF3)ppy

        Args:
            ligand_str: String containing ligand specifications

        Returns:
            List of ParsedLigand objects
        """
        ligands: List[ParsedLigand] = []
        i = 0

        while i < len(ligand_str):
            # Skip whitespace
            if ligand_str[i].isspace():
                i += 1
                continue

            # Handle parenthesized ligands
            if ligand_str[i] == "(":
                ligand_name, count, new_i = self._parse_parenthesized_ligand(
                    ligand_str, i
                )
                ligands.append(ParsedLigand(name=ligand_name, count=count))
                i = new_i
            else:
                # Try to match a known ligand
                matched, new_i = self._try_match_known_ligand(ligand_str, i, ligands)
                if matched:
                    i = new_i
                else:
                    # Extract unknown token
                    token, new_i = self._extract_unknown_token(ligand_str, i)
                    if token:
                        name, count = self._split_name_and_count(token)
                        ligands.append(ParsedLigand(name=name, count=count))
                    i = new_i

        return ligands

    def _parse_parenthesized_ligand(
        self, ligand_str: str, start: int
    ) -> Tuple[str, int, int]:
        """
        Parse a parenthesized ligand expression.

        Args:
            ligand_str: Full ligand string
            start: Index of opening parenthesis

        Returns:
            Tuple of (ligand_name, count, new_position)
        """
        # Find matching closing parenthesis
        depth = 1
        j = start + 1
        while j < len(ligand_str) and depth > 0:
            if ligand_str[j] == "(":
                depth += 1
            elif ligand_str[j] == ")":
                depth -= 1
            j += 1

        ligand_name = ligand_str[start + 1 : j - 1]

        # Check for count after closing parenthesis
        count = 1
        if j < len(ligand_str) and ligand_str[j].isdigit():
            count_str = ""
            while j < len(ligand_str) and ligand_str[j].isdigit():
                count_str += ligand_str[j]
                j += 1
            count = int(count_str)

        return ligand_name, count, j

    def _try_match_known_ligand(
        self, ligand_str: str, start: int, ligands: List[ParsedLigand]
    ) -> Tuple[bool, int]:
        """
        Try to match a known ligand at the current position.

        Args:
            ligand_str: Full ligand string
            start: Current position
            ligands: List to append matched ligand to

        Returns:
            Tuple of (matched: bool, new_position: int)
        """
        remaining = ligand_str[start:]

        # Sort by length (longest first) to match greedily
        for lig_name in sorted(self.ligand_db.keys(), key=len, reverse=True):
            if remaining.startswith(lig_name):
                j = start + len(lig_name)

                # Check for count after ligand
                count = 1
                if j < len(ligand_str) and ligand_str[j].isdigit():
                    count_str = ""
                    while j < len(ligand_str) and ligand_str[j].isdigit():
                        count_str += ligand_str[j]
                        j += 1
                    count = int(count_str)

                ligands.append(ParsedLigand(name=lig_name, count=count))
                return True, j

        return False, start

    def _extract_unknown_token(self, ligand_str: str, start: int) -> Tuple[str, int]:
        """
        Extract an unknown token until a delimiter is reached.

        Args:
            ligand_str: Full ligand string
            start: Current position

        Returns:
            Tuple of (token, new_position)
        """
        j = start
        delimiters = "()[]"

        while j < len(ligand_str) and ligand_str[j] not in delimiters:
            j += 1

        token = ligand_str[start:j] if j > start else ""
        return token, max(j, start + 1)

    def _split_name_and_count(self, token: str) -> Tuple[str, int]:
        """
        Split a token into name and trailing count.

        Args:
            token: Token string (e.g., "Cl2")

        Returns:
            Tuple of (name, count)
        """
        match = re.match(r"^(.+?)(\d+)$", token)
        if match:
            return match.group(1), int(match.group(2))
        return token, 1


# =============================================================================
# SMILES BUILDER
# =============================================================================


class SMILESBuilderError(Exception):
    """Custom exception for SMILES building errors."""

    pass


class SMILESBuilder:
    """
    Builds SMILES strings from parsed complex data.

    Note:
        SMILES has inherent limitations for representing coordination
        compounds. This builder produces a valid SMILES that captures
        the molecular components but may not represent the true
        bonding topology or stereochemistry.

    Example:
        >>> builder = SMILESBuilder()
        >>> parsed = ParsedComplex(metal="Ir", ligands=[...], ...)
        >>> smiles = builder.build(parsed)
    """

    def __init__(
        self,
        ligand_db: Optional[Dict[str, LigandInfo]] = None,
        metal_db: Optional[Dict[str, MetalInfo]] = None,
        counter_ion_db: Optional[Dict[str, LigandInfo]] = None,
    ) -> None:
        """
        Initialize the SMILES builder with chemical databases.

        Args:
            ligand_db: Ligand database. Uses LIGAND_DATABASE if None.
            metal_db: Metal database. Uses METAL_DATABASE if None.
            counter_ion_db: Counter ion database. Uses COUNTER_ION_DATABASE if None.
        """
        self.ligand_db = ligand_db if ligand_db is not None else LIGAND_DATABASE
        self.metal_db = metal_db if metal_db is not None else METAL_DATABASE
        self.counter_ion_db = (
            counter_ion_db if counter_ion_db is not None else COUNTER_ION_DATABASE
        )

    def build(self, parsed: ParsedComplex) -> str:
        """
        Build a SMILES string from a parsed complex.

        Args:
            parsed: ParsedComplex object

        Returns:
            SMILES string representation

        Raises:
            SMILESBuilderError: If SMILES cannot be constructed
        """
        # Calculate metal oxidation state
        metal_charge = self._calculate_metal_charge(parsed)

        # Build metal center SMILES
        metal_smiles = self._format_metal_smiles(parsed.metal, metal_charge)

        # Build ligand SMILES
        ligand_smiles_list = self._build_ligand_smiles(parsed.ligands)

        # Combine into single complex unit
        complex_parts = [metal_smiles] + ligand_smiles_list
        single_unit_smiles = ".".join(complex_parts)

        # Handle multiplicity (dimers, etc.)
        if parsed.multiplicity > 1:
            full_smiles = ".".join([single_unit_smiles] * parsed.multiplicity)
        else:
            full_smiles = single_unit_smiles

        # Add counter ions
        for ion_name, ion_count in parsed.counter_ions:
            ion_smiles = self._get_counter_ion_smiles(ion_name)
            for _ in range(ion_count):
                full_smiles += "." + ion_smiles

        return full_smiles

    def _calculate_metal_charge(self, parsed: ParsedComplex) -> int:
        """
        Calculate the formal oxidation state of the metal center.

        Uses charge balance:
            metal_charge = complex_charge - sum(ligand_charges)

        Args:
            parsed: ParsedComplex object

        Returns:
            Calculated metal oxidation state
        """
        # Sum up ligand charges
        total_ligand_charge = 0
        for ligand in parsed.ligands:
            if ligand.name in self.ligand_db:
                lig_info = self.ligand_db[ligand.name]
                total_ligand_charge += lig_info.charge * ligand.count

        # Determine complex charge
        if parsed.counter_ions:
            # If counter ions present, calculate complex charge from them
            counter_ion_charge = 0
            for ion_name, ion_count in parsed.counter_ions:
                if ion_name in self.counter_ion_db:
                    counter_ion_charge += (
                        self.counter_ion_db[ion_name].charge * ion_count
                    )
            # Overall compound is neutral: complex_charge + counter_ion_charge = 0
            complex_charge = -counter_ion_charge
        else:
            complex_charge = parsed.complex_charge

        # Adjust for multiplicity
        complex_charge_per_unit = complex_charge // parsed.multiplicity

        # Calculate metal charge
        metal_charge = complex_charge_per_unit - total_ligand_charge

        return metal_charge

    def _format_metal_smiles(self, symbol: str, charge: int) -> str:
        """
        Format a metal symbol with charge as SMILES.

        Args:
            symbol: Metal element symbol
            charge: Formal charge

        Returns:
            SMILES representation (e.g., "[Ir+3]")
        """
        if charge == 0:
            return f"[{symbol}]"
        elif charge > 0:
            if charge == 1:
                return f"[{symbol}+]"
            else:
                return f"[{symbol}+{charge}]"
        else:  # charge < 0
            if charge == -1:
                return f"[{symbol}-]"
            else:
                return f"[{symbol}{charge}]"

    def _build_ligand_smiles(self, ligands: List[ParsedLigand]) -> List[str]:
        """
        Build SMILES strings for all ligands.

        Args:
            ligands: List of ParsedLigand objects

        Returns:
            List of SMILES strings
        """
        smiles_list: List[str] = []

        for ligand in ligands:
            lig_smiles = self._get_ligand_smiles(ligand.name)
            # Add SMILES for each copy of the ligand
            for _ in range(ligand.count):
                smiles_list.append(lig_smiles)

        return smiles_list

    def _get_ligand_smiles(self, name: str) -> str:
        """
        Get SMILES for a ligand by name.

        Attempts:
            1. Direct lookup in database
            2. Search aliases
            3. Match as suffix (for modified ligands)

        Args:
            name: Ligand name/abbreviation

        Returns:
            SMILES string

        Raises:
            SMILESBuilderError: If ligand not found
        """
        # Direct lookup
        if name in self.ligand_db:
            return self.ligand_db[name].smiles

        # Check aliases
        for lig_name, info in self.ligand_db.items():
            if name in info.aliases:
                return info.smiles

        # Try to match as a modified ligand (e.g., "dF(CF3)ppy" -> use "ppy")
        for base_name in sorted(self.ligand_db.keys(), key=len, reverse=True):
            if name.endswith(base_name):
                # Found base ligand with modifications
                # TODO: In future, handle modifications properly
                return self.ligand_db[base_name].smiles

        raise SMILESBuilderError(f"Unknown ligand: '{name}'")

    def _get_counter_ion_smiles(self, name: str) -> str:
        """
        Get SMILES for a counter ion by name.

        Args:
            name: Counter ion name

        Returns:
            SMILES string

        Raises:
            SMILESBuilderError: If counter ion not found
        """
        if name in self.counter_ion_db:
            return self.counter_ion_db[name].smiles

        # Check aliases
        for ion_name, info in self.counter_ion_db.items():
            if name in info.aliases:
                return info.smiles

        raise SMILESBuilderError(f"Unknown counter ion: '{name}'")


# =============================================================================
# MAIN CONVERTER CLASS
# =============================================================================


class InorganicNameToSMILES:
    """
    Main converter class for inorganic/organometallic names to SMILES.

    This class provides a high-level interface for converting coordination
    complex names to SMILES notation.

    Example:
        >>> converter = InorganicNameToSMILES()
        >>> smiles = converter.convert("[IrCl(cod)]2")
        >>> print(smiles)
        '[Ir+].[Cl-].C1=CCCC=CCC1.[Ir+].[Cl-].C1=CCCC=CCC1'
    """

    def __init__(
        self,
        ligand_db: Optional[Dict[str, LigandInfo]] = None,
        metal_db: Optional[Dict[str, MetalInfo]] = None,
        counter_ion_db: Optional[Dict[str, LigandInfo]] = None,
    ) -> None:
        """
        Initialize the converter with optional custom databases.

        Args:
            ligand_db: Custom ligand database (optional)
            metal_db: Custom metal database (optional)
            counter_ion_db: Custom counter ion database (optional)
        """
        self.ligand_db = ligand_db if ligand_db is not None else LIGAND_DATABASE
        self.metal_db = metal_db if metal_db is not None else METAL_DATABASE
        self.counter_ion_db = (
            counter_ion_db if counter_ion_db is not None else COUNTER_ION_DATABASE
        )

        self.parser = ComplexNameParser(
            self.ligand_db, self.metal_db, self.counter_ion_db
        )
        self.builder = SMILESBuilder(self.ligand_db, self.metal_db, self.counter_ion_db)

    def convert(self, name: str) -> str:
        """
        Convert an inorganic complex name to SMILES.

        Args:
            name: The complex name (e.g., "[IrCl(cod)]2")

        Returns:
            SMILES string representation

        Raises:
            ParserError: If the name cannot be parsed
            SMILESBuilderError: If SMILES cannot be constructed
        """
        parsed = self.parser.parse(name)
        return self.builder.build(parsed)

    def convert_with_details(self, name: str) -> Tuple[str, ParsedComplex]:
        """
        Convert and return both SMILES and parsed structure.

        Args:
            name: The complex name

        Returns:
            Tuple of (SMILES string, ParsedComplex object)
        """
        parsed = self.parser.parse(name)
        smiles = self.builder.build(parsed)
        return smiles, parsed

    def add_ligand(
        self,
        name: str,
        smiles: str,
        denticity: int = 1,
        charge: int = 0,
        aliases: Optional[Tuple[str, ...]] = None,
        description: str = "",
    ) -> None:
        """
        Add a new ligand to the database.

        Args:
            name: Ligand abbreviation (e.g., "dppe")
            smiles: SMILES representation
            denticity: Number of coordination sites
            charge: Formal charge of the ligand
            aliases: Alternative names
            description: Human-readable description
        """
        self.ligand_db[name] = LigandInfo(
            smiles=smiles,
            denticity=denticity,
            charge=charge,
            aliases=aliases if aliases else tuple(),
            description=description,
        )

    def add_counter_ion(
        self,
        name: str,
        smiles: str,
        charge: int = -1,
        aliases: Optional[Tuple[str, ...]] = None,
        description: str = "",
    ) -> None:
        """
        Add a new counter ion to the database.

        Args:
            name: Counter ion abbreviation
            smiles: SMILES representation
            charge: Formal charge (usually -1)
            aliases: Alternative names
            description: Human-readable description
        """
        self.counter_ion_db[name] = LigandInfo(
            smiles=smiles,
            charge=charge,
            aliases=aliases if aliases else tuple(),
            description=description,
        )

    def list_available_ligands(self) -> List[str]:
        """Return list of available ligand abbreviations."""
        return list(self.ligand_db.keys())

    def list_available_metals(self) -> List[str]:
        """Return list of available metal symbols."""
        return list(self.metal_db.keys())

    def list_available_counter_ions(self) -> List[str]:
        """Return list of available counter ion abbreviations."""
        return list(self.counter_ion_db.keys())


# =============================================================================
# TESTING / DEMONSTRATION
# =============================================================================


def run_tests() -> None:
    """Run demonstration tests for the converter."""
    converter = InorganicNameToSMILES()

    test_cases = [
        "[IrCl(cod)]2",
        "[RhCp*Cl2]2",
        "[Ir(ppy)2(bpy)]PF6",
        "[Pd(PPh3)4]",
        "[PtCl2(en)]",
        "[Ru(bpy)3]2+",
        "[Fe(CO)5]",
    ]

    print("=" * 70)
    print("Inorganic Name to SMILES Converter - Test Results")
    print("=" * 70)

    for name in test_cases:
        print(f"\n{'─' * 70}")
        print(f"Input: {name}")
        print("─" * 70)

        try:
            smiles, parsed = converter.convert_with_details(name)

            print(f"  Metal:        {parsed.metal}")
            print(
                f"  Ligands:      {[(lig.name, lig.count) for lig in parsed.ligands]}"
            )
            print(f"  Multiplicity: {parsed.multiplicity}")
            print(f"  Charge:       {parsed.complex_charge}")
            print(f"  Counter ions: {parsed.counter_ions}")
            print(f"  SMILES:       {smiles}")

        except (ParserError, SMILESBuilderError) as e:
            print(f"  ERROR: {type(e).__name__}: {e}")

    print("\n" + "=" * 70)
    print("Available ligands:", converter.list_available_ligands())
    print("Available metals:", converter.list_available_metals())
    print("=" * 70)


if __name__ == "__main__":
    run_tests()
