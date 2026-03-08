from dataclasses import dataclass, field
from enum import Enum, auto
from typing import Dict, Tuple, TypedDict


class LigandType(Enum):
    """Classification of ligand charge types."""

    NEUTRAL = auto()
    ANIONIC = auto()
    CATIONIC = auto()


class BondingMode(TypedDict):
    name: str
    coordination_kind: str
    hapticity: int | None
    donor_mapnums: Tuple[int, ...]
    preferred_bond_type: str


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
    mapped_smiles: str
    denticity: int = 1
    charge: int = 0
    aliases: Tuple[str, ...] = field(default_factory=tuple)
    description: str = ""
    binding_modes: Tuple[BondingMode, ...] = field(default_factory=tuple)

    @property
    def ligand_type(self) -> LigandType:
        """Determine ligand type based on formal charge."""
        if self.charge < 0:
            return LigandType.ANIONIC
        elif self.charge > 0:
            return LigandType.CATIONIC
        return LigandType.NEUTRAL

    @property
    def rdkit_charge(self) -> int:
        """Get the RDKit charge of the ligand."""
        from rdkit import Chem

        mol = Chem.MolFromSmiles(self.smiles)
        if mol is None:
            return 0
        return Chem.GetFormalCharge(mol)


@dataclass
class MetalInfo:
    """
    Information about a transition metal.

    Attributes:
        symbol: Element symbol (e.g., "Ir")
        name: Full element name (e.g., "Iridium")
        common_oxidation_states: List of typical oxidation states
        atomic_number: Atomic number of the element
        bond_type: Type of bond (e.g., "single", "double", "triple")
        bond_atom_idx: Index of atoms involved in the bond
    """

    symbol: str
    name: str
    common_oxidation_states: Tuple[int, ...]
    atomic_number: int


LIGAND_DATABASE: Dict[str, LigandInfo] = {
    # -------------------------------------------------------------------------
    # Monodentate Neutral Ligands
    # -------------------------------------------------------------------------
    "CO": LigandInfo(
        smiles="[C-]#[O+]",
        mapped_smiles="[C-:1]#[O+:2]",
        denticity=1,
        charge=0,
        aliases=("carbonyl",),
        description="Carbonyl ligand",
        binding_modes=(
            {
                "name": "kappa-C",
                "coordination_kind": "kappa",
                "hapticity": None,
                "donor_mapnums": (1,),
                "preferred_bond_type": "dative",
            },
        ),
    ),
    "PPh3": LigandInfo(
        smiles="c1ccc(P(c2ccccc2)c2ccccc2)cc1",
        mapped_smiles="[cH:1]1[cH:2][cH:3][c:4]([P:5]([c:6]2[cH:7][cH:8][cH:9][cH:10][cH:11]2)[c:12]2[cH:13][cH:14][cH:15][cH:16][cH:17]2)[cH:18][cH:19]1",
        denticity=1,
        charge=0,
        aliases=("triphenylphosphine", "Ph3P"),
        description="Triphenylphosphine",
        binding_modes=(
            {
                "name": "kappa-P",
                "coordination_kind": "kappa",
                "hapticity": None,
                "donor_mapnums": (5,),
                "preferred_bond_type": "dative",
            },
        ),
    ),
    "py": LigandInfo(
        smiles="c1ccncc1",
        mapped_smiles="[cH:1]1[cH:2][cH:3][n:4][cH:5][cH:6]1",
        denticity=1,
        charge=0,
        aliases=("pyridine", "Py"),
        description="Pyridine",
        binding_modes=(
            {
                "name": "kappa-N",
                "coordination_kind": "kappa",
                "hapticity": None,
                "donor_mapnums": (4,),
                "preferred_bond_type": "dative",
            },
        ),
    ),
    "NH3": LigandInfo(
        smiles="N",
        mapped_smiles="[NH3:1]",
        denticity=1,
        charge=0,
        aliases=("ammonia", "ammine"),
        description="Ammonia/Ammine ligand",
        binding_modes=(
            {
                "name": "kappa-N",
                "coordination_kind": "kappa",
                "hapticity": None,
                "donor_mapnums": (1,),
                "preferred_bond_type": "dative",
            },
        ),
    ),
    "H2O": LigandInfo(
        smiles="O",
        mapped_smiles="[OH2:1]",
        denticity=1,
        charge=0,
        aliases=("water", "aqua", "aquo"),
        description="Water/Aqua ligand",
        binding_modes=(
            {
                "name": "kappa-O",
                "coordination_kind": "kappa",
                "hapticity": None,
                "donor_mapnums": (1,),
                "preferred_bond_type": "dative",
            },
        ),
    ),
    "MeCN": LigandInfo(
        smiles="CC#N",
        mapped_smiles="[CH3:1][C:2]#[N:3]",
        denticity=1,
        charge=0,
        aliases=("acetonitrile", "NCMe", "CH3CN"),
        description="Acetonitrile",
        binding_modes=(
            {
                "name": "kappa-N",
                "coordination_kind": "kappa",
                "hapticity": None,
                "donor_mapnums": (3,),
                "preferred_bond_type": "dative",
            },
        ),
    ),
    "PMe3": LigandInfo(
        smiles="",
        mapped_smiles="",
        denticity=1,
        charge=0,
        aliases=("trimethylphosphine",),
        description="Trimethylphosphine",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "PEt3": LigandInfo(
        smiles="CCP(CC)CC",
        mapped_smiles="[CH3:1][CH2:2][P:3]([CH2:4][CH3:5])[CH2:6][CH3:7]",
        denticity=1,
        charge=0,
        aliases=("triethylphosphine",),
        description="Triethylphosphine",
        binding_modes=(
            {
                "name": "kappa-P",
                "coordination_kind": "kappa",
                "hapticity": None,
                "donor_mapnums": (3,),
                "preferred_bond_type": "dative",
            },
        ),
    ),
    "PiPr3": LigandInfo(
        smiles="CC(C)P(C(C)C)C(C)C",
        mapped_smiles="[CH3:1][CH:2]([CH3:3])[P:4]([CH:5]([CH3:6])[CH3:7])[CH:8]([CH3:9])[CH3:10]",
        denticity=1,
        charge=0,
        aliases=("triisopropylphosphine", "P(iPr)3"),
        description="Triisopropylphosphine",
        binding_modes=(
            {
                "name": "kappa-P",
                "coordination_kind": "kappa",
                "hapticity": None,
                "donor_mapnums": (4,),
                "preferred_bond_type": "dative",
            },
        ),
    ),
    "PCy3": LigandInfo(
        smiles="C1CCC(P(C2CCCCC2)C2CCCCC2)CC1",
        mapped_smiles="[CH2:1]1[CH2:2][CH2:3][CH:4]([P:5]([CH:6]2[CH2:7][CH2:8][CH2:9][CH2:10][CH2:11]2)[CH:12]2[CH2:13][CH2:14][CH2:15][CH2:16][CH2:17]2)[CH2:18][CH2:19]1",
        denticity=1,
        charge=0,
        aliases=("tricyclohexylphosphine",),
        description="Tricyclohexylphosphine",
        binding_modes=(
            {
                "name": "kappa-P",
                "coordination_kind": "kappa",
                "hapticity": None,
                "donor_mapnums": (5,),
                "preferred_bond_type": "dative",
            },
        ),
    ),
    "PtBu3": LigandInfo(
        smiles="CCCCP(CCCC)CCCC",
        mapped_smiles="[CH3:1][CH2:2][CH2:3][CH2:4][P:5]([CH2:6][CH2:7][CH2:8][CH3:9])[CH2:10][CH2:11][CH2:12][CH3:13]",
        denticity=1,
        charge=0,
        aliases=("tri-tert-butylphosphine", "P(tBu)3"),
        description="Tri-tert-butylphosphine",
        binding_modes=(
            {
                "name": "kappa-P",
                "coordination_kind": "kappa",
                "hapticity": None,
                "donor_mapnums": (5,),
                "preferred_bond_type": "dative",
            },
        ),
    ),
    "P(OMe)3": LigandInfo(
        smiles="CO[PH](C)(OC)OC",
        mapped_smiles="[CH3:1][O:2][PH:3]([CH3:4])([O:5][CH3:6])[O:7][CH3:8]",
        denticity=1,
        charge=0,
        aliases=("trimethylphosphite",),
        description="Trimethyl phosphite",
        binding_modes=(
            {
                "name": "kappa-P",
                "coordination_kind": "kappa",
                "hapticity": None,
                "donor_mapnums": (3,),
                "preferred_bond_type": "dative",
            },
        ),
    ),
    "P(OEt)3": LigandInfo(
        smiles="CCOP(OCC)OCC",
        mapped_smiles="[CH3:1][CH2:2][O:3][P:4]([O:5][CH2:6][CH3:7])[O:8][CH2:9][CH3:10]",
        denticity=1,
        charge=0,
        aliases=("triethylphosphite",),
        description="Triethyl phosphite",
        binding_modes=(
            {
                "name": "kappa-P",
                "coordination_kind": "kappa",
                "hapticity": None,
                "donor_mapnums": (4,),
                "preferred_bond_type": "dative",
            },
        ),
    ),
    "P(OPh)3": LigandInfo(
        smiles="c1ccc(OP(Oc2ccccc2)Oc2ccccc2)cc1",
        mapped_smiles="[cH:1]1[cH:2][cH:3][c:4]([O:5][P:6]([O:7][c:8]2[cH:9][cH:10][cH:11][cH:12][cH:13]2)[O:14][c:15]2[cH:16][cH:17][cH:18][cH:19][cH:20]2)[cH:21][cH:22]1",
        denticity=1,
        charge=0,
        aliases=("triphenylphosphite",),
        description="Triphenyl phosphite",
        binding_modes=(
            {
                "name": "kappa-P",
                "coordination_kind": "kappa",
                "hapticity": None,
                "donor_mapnums": (6,),
                "preferred_bond_type": "dative",
            },
        ),
    ),
    # N-Heterocyclic Carbenes (NHCs) - hugely important in modern catalysis
    "IMes": LigandInfo(
        smiles="Cc1cc(C)c(N2[C]N(c3c(C)cc(C)cc3C)C=C2)c(C)c1",
        mapped_smiles="[CH3:1][c:2]1[cH:3][c:4]([CH3:5])[c:6]([N:7]2[C:8][N:9]([c:10]3[c:11]([CH3:12])[cH:13][c:14]([CH3:15])[cH:16][c:17]3[CH3:18])[CH:19]=[CH:20]2)[c:21]([CH3:22])[cH:23]1",
        denticity=1,
        charge=0,
        aliases=("1,3-bis(2,4,6-trimethylphenyl)imidazol-2-ylidene",),
        description="IMes carbene",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "IPr": LigandInfo(
        smiles="CC(C)c1cccc(C(C)C)c1N1[C]N(c2c(C(C)C)cccc2C(C)C)C=C1",
        mapped_smiles="[CH3:1][CH:2]([CH3:3])[c:4]1[cH:5][cH:6][cH:7][c:8]([CH:9]([CH3:10])[CH3:11])[c:12]1[N:13]1[C:14][N:15]([c:16]2[c:17]([CH:18]([CH3:19])[CH3:20])[cH:21][cH:22][cH:23][c:24]2[CH:25]([CH3:26])[CH3:27])[CH:28]=[CH:29]1",
        denticity=1,
        charge=0,
        aliases=("1,3-bis(2,6-diisopropylphenyl)imidazol-2-ylidene",),
        description="IPr carbene",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "SIMes": LigandInfo(
        smiles="Cc1cc(C)c(N2[C]N(c3c(C)cc(C)cc3C)CC2)c(C)c1",
        mapped_smiles="[CH3:1][c:2]1[cH:3][c:4]([CH3:5])[c:6]([N:7]2[C:8][N:9]([c:10]3[c:11]([CH3:12])[cH:13][c:14]([CH3:15])[cH:16][c:17]3[CH3:18])[CH2:19][CH2:20]2)[c:21]([CH3:22])[cH:23]1",
        denticity=1,
        charge=0,
        aliases=("1,3-bis(2,4,6-trimethylphenyl)imidazolidin-2-ylidene",),
        description="Saturated IMes carbene",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "SIPr": LigandInfo(
        smiles="CC(C)c1cccc(C(C)C)c1N1[C-]=[N+](c2c(C(C)C)cccc2C(C)C)CC1",
        mapped_smiles="[CH3:1][CH:2]([CH3:3])[c:4]1[cH:5][cH:6][cH:7][c:8]([CH:9]([CH3:10])[CH3:11])[c:12]1[N:13]1[C-:14]=[N+:15]([c:16]2[c:17]([CH:18]([CH3:19])[CH3:20])[cH:21][cH:22][cH:23][c:24]2[CH:25]([CH3:26])[CH3:27])[CH2:28][CH2:29]1",
        denticity=1,
        charge=0,
        aliases=("1,3-bis(2,6-diisopropylphenyl)imidazolidin-2-ylidene",),
        description="Saturated IPr carbene",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "ICy": LigandInfo(
        smiles="[C]1N(C2CCCCC2)C=CN1C1CCCCC1",
        mapped_smiles="[C:1]1[N:2]([CH:3]2[CH2:4][CH2:5][CH2:6][CH2:7][CH2:8]2)[CH:9]=[CH:10][N:11]1[CH:12]1[CH2:13][CH2:14][CH2:15][CH2:16][CH2:17]1",
        denticity=1,
        charge=0,
        aliases=("1,3-dicyclohexylimidazol-2-ylidene",),
        description="ICy carbene",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "ItBu": LigandInfo(
        smiles="CC(C)(C)N1[C]N(C(C)(C)C)C=C1",
        mapped_smiles="[CH3:1][C:2]([CH3:3])([CH3:4])[N:5]1[C:6][N:7]([C:8]([CH3:9])([CH3:10])[CH3:11])[CH:12]=[CH:13]1",
        denticity=1,
        charge=0,
        aliases=("1,3-di-tert-butylimidazol-2-ylidene",),
        description="ItBu carbene",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "IMe": LigandInfo(
        smiles="CN1[C]N(C)C=C1",
        mapped_smiles="[CH3:1][N:2]1[C:3][N:4]([CH3:5])[CH:6]=[CH:7]1",
        denticity=1,
        charge=0,
        aliases=("1,3-dimethylimidazol-2-ylidene",),
        description="IMe carbene",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "IAd": LigandInfo(
        smiles="[C]1N(C23CC4CC(CC(C4)C2)C3)C=CN1C12CC3CC(CC(C3)C1)C2",
        mapped_smiles="[C:1]1[N:2]([C:3]23[CH2:4][CH:5]4[CH2:6][CH:7]([CH2:8][CH:9]([CH2:10]4)[CH2:11]2)[CH2:12]3)[CH:13]=[CH:14][N:15]1[C:16]12[CH2:17][CH:18]3[CH2:19][CH:20]([CH2:21][CH:22]([CH2:23]3)[CH2:24]1)[CH2:25]2",
        denticity=1,
        charge=0,
        aliases=("1,3-bis(adamantyl)imidazol-2-ylidene",),
        description="IAd carbene",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    # Other common neutral donors
    "THF": LigandInfo(
        smiles="C1CCOC1",
        mapped_smiles="[CH2:1]1[CH2:2][CH2:3][O:4][CH2:5]1",
        denticity=1,
        charge=0,
        aliases=("tetrahydrofuran", "thf"),
        description="Tetrahydrofuran",
        binding_modes=(
            {
                "name": "kappa-O",
                "coordination_kind": "kappa",
                "hapticity": None,
                "donor_mapnums": (4,),
                "preferred_bond_type": "dative",
            },
        ),
    ),
    "Et2O": LigandInfo(
        smiles="CCOCC",
        mapped_smiles="[CH3:1][CH2:2][O:3][CH2:4][CH3:5]",
        denticity=1,
        charge=0,
        aliases=("diethylether", "ether"),
        description="Diethyl ether",
        binding_modes=(
            {
                "name": "kappa-O",
                "coordination_kind": "kappa",
                "hapticity": None,
                "donor_mapnums": (3,),
                "preferred_bond_type": "dative",
            },
        ),
    ),
    "DMF": LigandInfo(
        smiles="CN(C)C=O",
        mapped_smiles="[CH3:1][N:2]([CH3:3])[CH:4]=[O:5]",
        denticity=1,
        charge=0,
        aliases=("dimethylformamide",),
        description="Dimethylformamide",
        binding_modes=(
            {
                "name": "kappa-N",
                "coordination_kind": "kappa",
                "hapticity": None,
                "donor_mapnums": (2,),
                "preferred_bond_type": "dative",
            },
        ),
    ),
    "DMSO": LigandInfo(
        smiles="CS(C)=O",
        mapped_smiles="[CH3:1][S:2]([CH3:3])=[O:4]",
        denticity=1,
        charge=0,
        aliases=("dimethylsulfoxide",),
        description="Dimethyl sulfoxide",
        binding_modes=(
            {
                "name": "kappa-S",
                "coordination_kind": "kappa",
                "hapticity": None,
                "donor_mapnums": (2,),
                "preferred_bond_type": "dative",
            },
        ),
    ),
    "NMe3": LigandInfo(
        smiles="CN(C)C",
        mapped_smiles="[CH3:1][N:2]([CH3:3])[CH3:4]",
        denticity=1,
        charge=0,
        aliases=("trimethylamine",),
        description="Trimethylamine",
        binding_modes=(
            {
                "name": "kappa-N",
                "coordination_kind": "kappa",
                "hapticity": None,
                "donor_mapnums": (2,),
                "preferred_bond_type": "dative",
            },
        ),
    ),
    "NEt3": LigandInfo(
        smiles="CCN(CC)CC",
        mapped_smiles="[CH3:1][CH2:2][N:3]([CH2:4][CH3:5])[CH2:6][CH3:7]",
        denticity=1,
        charge=0,
        aliases=("triethylamine", "TEA"),
        description="Triethylamine",
        binding_modes=(
            {
                "name": "kappa-N",
                "coordination_kind": "kappa",
                "hapticity": None,
                "donor_mapnums": (3,),
                "preferred_bond_type": "dative",
            },
        ),
    ),
    "DMAP": LigandInfo(
        smiles="Cc1ccccc1N(C)C",
        mapped_smiles="[CH3:1][c:2]1[cH:3][cH:4][cH:5][cH:6][c:7]1[N:8]([CH3:9])[CH3:10]",
        denticity=1,
        charge=0,
        aliases=("4-dimethylaminopyridine",),
        description="4-Dimethylaminopyridine",
        binding_modes=(
            {
                "name": "kappa-N",
                "coordination_kind": "kappa",
                "hapticity": None,
                "donor_mapnums": (8,),
                "preferred_bond_type": "dative",
            },
        ),
    ),
    # Isocyanides
    "CNtBu": LigandInfo(
        smiles="CC(C)(C)C#N",
        mapped_smiles="[CH3:1][C:2]([CH3:3])([CH3:4])[C:5]#[N:6]",
        denticity=1,
        charge=0,
        aliases=("tert-butylisocyanide",),
        description="tert-Butyl isocyanide",
        binding_modes=(
            {
                "name": "kappa-C",
                "coordination_kind": "kappa",
                "hapticity": None,
                "donor_mapnums": (5,),
                "preferred_bond_type": "dative",
            },
        ),
    ),
    "CNXyl": LigandInfo(
        smiles="[C-]#[N+]c1c(C)cccc1C",
        mapped_smiles="[C-:1]#[N+:2][c:3]1[c:4]([CH3:5])[cH:6][cH:7][cH:8][c:9]1[CH3:10]",
        denticity=1,
        charge=0,
        aliases=("2,6-xylylisocyanide", "2,6-dimethylphenylisocyanide"),
        description="2,6-Xylyl isocyanide",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "CNPh": LigandInfo(
        smiles="O=C=Nc1ccccc1",
        mapped_smiles="[O:1]=[C:2]=[N:3][c:4]1[cH:5][cH:6][cH:7][cH:8][cH:9]1",
        denticity=1,
        charge=0,
        aliases=("phenylisocyanide",),
        description="Phenyl isocyanide",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "CNMe": LigandInfo(
        smiles="[C-]#[N+]C",
        mapped_smiles="[C-:1]#[N+:2][CH3:3]",
        denticity=1,
        charge=0,
        aliases=("methylisocyanide",),
        description="Methyl isocyanide",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    # Carbene ligands (non-NHC)
    "CHPh": LigandInfo(
        smiles="[CH]c1ccccc1",
        mapped_smiles="[CH:1][c:2]1[cH:3][cH:4][cH:5][cH:6][cH:7]1",
        denticity=1,
        charge=0,
        aliases=("benzylidene", "phenylcarbene"),
        description="Benzylidene",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    # Olefins (η²)
    "ethylene": LigandInfo(
        smiles="C=C",
        mapped_smiles="[CH2:1]=[CH2:2]",
        denticity=1,
        charge=0,
        aliases=("C2H4", "eth"),
        description="Ethylene (η²)",
        binding_modes=(
            {
                "name": "eta2",
                "coordination_kind": "haptic",
                "hapticity": 2,
                "donor_mapnums": (1, 2),
                "preferred_bond_type": "dative",
            },
        ),
    ),
    # Dinitrogen and other small molecules
    "N2": LigandInfo(
        smiles="N#N",
        mapped_smiles="[N:1]#[N:2]",
        denticity=1,
        charge=0,
        aliases=("dinitrogen",),
        description="Dinitrogen",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "dative",
            },
        ),
    ),
    "NO": LigandInfo(
        smiles="N=O",
        mapped_smiles="[NH:1]=[O:2]",
        denticity=1,
        charge=0,
        aliases=("nitrosyl",),
        description="Nitrosyl (neutral counting)",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "dative",
            },
        ),
    ),
    "CS": LigandInfo(
        smiles="C=S",
        mapped_smiles="[CH2:1]=[S:2]",
        denticity=1,
        charge=0,
        aliases=("thiocarbonyl",),
        description="Thiocarbonyl",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "dative",
            },
        ),
    ),
    "SO2": LigandInfo(
        smiles="O=S=O",
        mapped_smiles="[O:1]=[S:2]=[O:3]",
        denticity=1,
        charge=0,
        aliases=("sulfurdioxide",),
        description="Sulfur dioxide",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "dative",
            },
        ),
    ),
    "O2": LigandInfo(
        smiles="O=O",
        mapped_smiles="[O:1]=[O:2]",
        denticity=1,
        charge=0,
        aliases=("dioxygen",),
        description="Dioxygen",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "dative",
            },
        ),
    ),
    "H2": LigandInfo(
        smiles="",
        mapped_smiles="",
        denticity=1,
        charge=0,
        aliases=("dihydrogen",),
        description="Dihydrogen (η²-H2)",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "dative",
            },
        ),
    ),
    # -------------------------------------------------------------------------
    # Monodentate Anionic Ligands
    # -------------------------------------------------------------------------
    "Cl": LigandInfo(
        smiles="[Cl-]",
        mapped_smiles="[Cl-:1]",
        denticity=1,
        charge=-1,
        aliases=("chloro", "chloride", "chlorido"),
        description="Chloride",
        binding_modes=(
            {
                "name": "kappa-Cl",
                "coordination_kind": "kappa",
                "hapticity": None,
                "donor_mapnums": (1,),
                "preferred_bond_type": "single",
            },
        ),
    ),
    "Br": LigandInfo(
        smiles="[Br-]",
        mapped_smiles="[Br-:1]",
        denticity=1,
        charge=-1,
        aliases=("bromo", "bromide", "bromido"),
        description="Bromide",
        binding_modes=(
            {
                "name": "kappa-Br",
                "coordination_kind": "kappa",
                "hapticity": None,
                "donor_mapnums": (1,),
                "preferred_bond_type": "single",
            },
        ),
    ),
    "I": LigandInfo(
        smiles="[I-]",
        mapped_smiles="[I-:1]",
        denticity=1,
        charge=-1,
        aliases=("iodo", "iodide", "iodido"),
        description="Iodide",
        binding_modes=(
            {
                "name": "kappa-I",
                "coordination_kind": "kappa",
                "hapticity": None,
                "donor_mapnums": (1,),
                "preferred_bond_type": "single",
            },
        ),
    ),
    "F": LigandInfo(
        smiles="[F-]",
        mapped_smiles="[F-:1]",
        denticity=1,
        charge=-1,
        aliases=("fluoro", "fluoride", "fluorido"),
        description="Fluoride",
        binding_modes=(
            {
                "name": "kappa-F",
                "coordination_kind": "kappa",
                "hapticity": None,
                "donor_mapnums": (1,),
                "preferred_bond_type": "single",
            },
        ),
    ),
    "H": LigandInfo(
        smiles="[H-]",
        mapped_smiles="[H-:1]",
        denticity=1,
        charge=-1,
        aliases=("hydrido", "hydride"),
        description="Hydride",
        binding_modes=(
            {
                "name": "kappa-H",
                "coordination_kind": "kappa",
                "hapticity": None,
                "donor_mapnums": (1,),
                "preferred_bond_type": "single",
            },
        ),
    ),
    "CN": LigandInfo(
        smiles="[C-]#N",
        mapped_smiles="[C-:1]#[N:2]",
        denticity=1,
        charge=-1,
        aliases=("cyano", "cyanide", "cyanido"),
        description="Cyanide",
        binding_modes=(
            {
                "name": "kappa-C",
                "coordination_kind": "kappa",
                "hapticity": None,
                "donor_mapnums": (1,),
                "preferred_bond_type": "triple",
            },
        ),
    ),
    "OAc": LigandInfo(
        smiles="CC(=O)[O-]",
        mapped_smiles="[CH3:1][C:2](=[O:3])[O-:4]",
        denticity=1,
        charge=-1,
        aliases=("acetato", "acetate", "OAc"),
        description="Acetate",
        binding_modes=(
            {
                "name": "kappa-O",
                "coordination_kind": "kappa",
                "hapticity": None,
                "donor_mapnums": (4,),
                "preferred_bond_type": "single",
            },
        ),
    ),
    "OMe": LigandInfo(
        smiles="C[O-]",
        mapped_smiles="[CH3:1][O-:2]",
        denticity=1,
        charge=-1,
        aliases=("methoxo", "methoxide"),
        description="Methoxide",
        binding_modes=(
            {
                "name": "kappa-O",
                "coordination_kind": "kappa",
                "hapticity": None,
                "donor_mapnums": (2,),
                "preferred_bond_type": "single",
            },
        ),
    ),
    # Common anions
    "OH": LigandInfo(
        smiles="[OH-]",
        mapped_smiles="[OH-:1]",
        denticity=1,
        charge=-1,
        aliases=("hydroxo", "hydroxide", "hydroxido"),
        description="Hydroxide",
        binding_modes=(
            {
                "name": "kappa-O",
                "coordination_kind": "kappa",
                "hapticity": None,
                "donor_mapnums": (1,),
                "preferred_bond_type": "single",
            },
        ),
    ),
    "OEt": LigandInfo(
        smiles="CC(C)[O-]",
        mapped_smiles="[CH3:1][CH:2]([CH3:3])[O-:4]",
        denticity=1,
        charge=-1,
        aliases=("ethoxo", "ethoxide"),
        description="Ethoxide",
        binding_modes=(
            {
                "name": "kappa-O",
                "coordination_kind": "kappa",
                "hapticity": None,
                "donor_mapnums": (4,),
                "preferred_bond_type": "single",
            },
        ),
    ),
    "OiPr": LigandInfo(
        smiles="CC(C)[O-]",
        mapped_smiles="[CH3:1][CH:2]([CH3:3])[O-:4]",
        denticity=1,
        charge=-1,
        aliases=("isopropoxo", "isopropoxide"),
        description="Isopropoxide",
        binding_modes=(
            {
                "name": "kappa-O",
                "coordination_kind": "kappa",
                "hapticity": None,
                "donor_mapnums": (4,),
                "preferred_bond_type": "single",
            },
        ),
    ),
    "OtBu": LigandInfo(
        smiles="CC(C)(C)[O-]",
        mapped_smiles="[CH3:1][C:2]([CH3:3])([CH3:4])[O-:5]",
        denticity=1,
        charge=-1,
        aliases=("tert-butoxo", "tert-butoxide"),
        description="tert-Butoxide",
        binding_modes=(
            {
                "name": "kappa-O",
                "coordination_kind": "kappa",
                "hapticity": None,
                "donor_mapnums": (5,),
                "preferred_bond_type": "single",
            },
        ),
    ),
    "OPh": LigandInfo(
        smiles="[O-]c1ccccc1",
        mapped_smiles="[O-:1][c:2]1[cH:3][cH:4][cH:5][cH:6][cH:7]1",
        denticity=1,
        charge=-1,
        aliases=("phenoxo", "phenoxide", "phenolate"),
        description="Phenoxide",
        binding_modes=(
            {
                "name": "kappa-O",
                "coordination_kind": "kappa",
                "hapticity": None,
                "donor_mapnums": (1,),
                "preferred_bond_type": "single",
            },
        ),
    ),
    # Alkyls
    "Me": LigandInfo(
        smiles="[CH3-]",
        mapped_smiles="[CH3-:1]",
        denticity=1,
        charge=-1,
        aliases=("methyl", "CH3"),
        description="Methyl",
        binding_modes=(
            {
                "name": "kappa-C",
                "coordination_kind": "kappa",
                "hapticity": None,
                "donor_mapnums": (1,),
                "preferred_bond_type": "single",
            },
        ),
    ),
    "Et": LigandInfo(
        smiles="[CH2-]C",
        mapped_smiles="[CH2-:1][CH3:2]",
        denticity=1,
        charge=-1,
        aliases=("ethyl", "C2H5"),
        description="Ethyl",
        binding_modes=(
            {
                "name": "kappa-C",
                "coordination_kind": "kappa",
                "hapticity": None,
                "donor_mapnums": (1,),
                "preferred_bond_type": "single",
            },
        ),
    ),
    # "nBu": LigandInfo(
    #     smiles="[CH2-]CC",
    #     mapped_smiles="[CH2-:1][CH2:2][CH3:3]",
    #     denticity=1,
    #     charge=-1,
    #     aliases=("n-butyl", "butyl"),
    #     description="n-Butyl",
    #     binding_modes=(
    #         {
    #             "name": "",
    #             "coordination_kind": "",
    #             "hapticity": None,
    #             "donor_mapnums": (),
    #             "preferred_bond_type": "",
    #         },
    #     ),
    # ),
    "Ph": LigandInfo(
        smiles="[c-]1ccccc1",
        mapped_smiles="[c-:1]1[cH:2][cH:3][cH:4][cH:5][cH:6]1",
        denticity=1,
        charge=-1,
        aliases=("phenyl", "C6H5"),
        description="Phenyl",
        binding_modes=(
            {
                "name": "kappa-C",
                "coordination_kind": "kappa",
                "hapticity": None,
                "donor_mapnums": (1,),
                "preferred_bond_type": "single",
            },
        ),
    ),
    "Bn": LigandInfo(
        smiles="[CH2-]c1ccccc1",
        mapped_smiles="[CH2-:1][c:2]1[cH:3][cH:4][cH:5][cH:6][cH:7]1",
        denticity=1,
        charge=-1,
        aliases=("benzyl",),
        description="Benzyl",
        binding_modes=(
            {
                "name": "kappa-C",
                "coordination_kind": "kappa",
                "hapticity": None,
                "donor_mapnums": (1,),
                "preferred_bond_type": "single",
            },
        ),
    ),
    "vinyl": LigandInfo(
        smiles="[C-]=C",
        mapped_smiles="[C-:1]=[CH2:2]",
        denticity=1,
        charge=-1,
        aliases=("ethenyl",),
        description="Vinyl",
        binding_modes=(
            {
                "name": "kappa-C",
                "coordination_kind": "kappa",
                "hapticity": None,
                "donor_mapnums": (1,),
                "preferred_bond_type": "single",
            },
        ),
    ),
    "allyl": LigandInfo(
        smiles="[CH-]C=C",
        mapped_smiles="[CH-:1][CH:2]=[CH2:3]",
        denticity=1,
        charge=-1,
        aliases=("η1-allyl", "propenyl", "C3H5"),
        description="Allyl (σ-bound)",
        binding_modes=(
            {
                "name": "kappa-C",
                "coordination_kind": "kappa",
                "hapticity": None,
                "donor_mapnums": (1,),
                "preferred_bond_type": "single",
            },
        ),
    ),
    "Np": LigandInfo(
        smiles="[CH2-]C(C)(C)C",
        mapped_smiles="[CH2-:1][C:2]([CH3:3])([CH3:4])[CH3:5]",
        denticity=1,
        charge=-1,
        aliases=("neopentyl",),
        description="Neopentyl",
        binding_modes=(
            {
                "name": "kappa-C",
                "coordination_kind": "kappa",
                "hapticity": None,
                "donor_mapnums": (1,),
                "preferred_bond_type": "single",
            },
        ),
    ),
    "Mes": LigandInfo(
        smiles="Cc1[c-]c(C)cc(C)c1",
        mapped_smiles="[CH3:1][c:2]1[c-:3][c:4]([CH3:5])[cH:6][c:7]([CH3:8])[cH:9]1",
        denticity=1,
        charge=-1,
        aliases=("mesityl", "2,4,6-trimethylphenyl"),
        description="Mesityl",
        binding_modes=(
            {
                "name": "kappa-C",
                "coordination_kind": "kappa",
                "hapticity": None,
                "donor_mapnums": (3,),
                "preferred_bond_type": "single",
            },
        ),
    ),
    # Silyls
    "SiMe3": LigandInfo(
        smiles="C[Si-](C)C",
        mapped_smiles="[CH3:1][Si-:2]([CH3:3])[CH3:4]",
        denticity=1,
        charge=-1,
        aliases=("trimethylsilyl", "TMS"),
        description="Trimethylsilyl",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "SiPh3": LigandInfo(
        smiles="c1ccc([Si-](c2ccccc2)c2ccccc2)cc1",
        mapped_smiles="[cH:1]1[cH:2][cH:3][c:4]([Si-:5]([c:6]2[cH:7][cH:8][cH:9][cH:10][cH:11]2)[c:12]2[cH:13][cH:14][cH:15][cH:16][cH:17]2)[cH:18][cH:19]1",
        denticity=1,
        charge=-1,
        aliases=("triphenylsilyl",),
        description="Triphenylsilyl",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    # Amides
    "NMe2": LigandInfo(
        smiles="C[N-]C",
        mapped_smiles="[CH3:1][N-:2][CH3:3]",
        denticity=1,
        charge=-1,
        aliases=("dimethylamido",),
        description="Dimethylamide",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "NEt2": LigandInfo(
        smiles="CC[N-]CC",
        mapped_smiles="[CH3:1][CH2:2][N-:3][CH2:4][CH3:5]",
        denticity=1,
        charge=-1,
        aliases=("diethylamido",),
        description="Diethylamide",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "NiPr2": LigandInfo(
        smiles="CC(C)[N-]C(C)C",
        mapped_smiles="[CH3:1][CH:2]([CH3:3])[N-:4][CH:5]([CH3:6])[CH3:7]",
        denticity=1,
        charge=-1,
        aliases=("diisopropylamido",),
        description="Diisopropylamide",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "NPh2": LigandInfo(
        smiles="c1ccc([N-]c2ccccc2)cc1",
        mapped_smiles="[cH:1]1[cH:2][cH:3][c:4]([N-:5][c:6]2[cH:7][cH:8][cH:9][cH:10][cH:11]2)[cH:12][cH:13]1",
        denticity=1,
        charge=-1,
        aliases=("diphenylamido",),
        description="Diphenylamide",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "NTMS2": LigandInfo(
        smiles="C[Si](C)(C)[N-][Si](C)(C)C",
        mapped_smiles="[CH3:1][Si:2]([CH3:3])([CH3:4])[N-:5][Si:6]([CH3:7])([CH3:8])[CH3:9]",
        denticity=1,
        charge=-1,
        aliases=("bis(trimethylsilyl)amido", "HMDS", "N(SiMe3)2"),
        description="Bis(trimethylsilyl)amide",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "NHPh": LigandInfo(
        smiles="[NH-]c1ccccc1",
        mapped_smiles="[NH-:1][c:2]1[cH:3][cH:4][cH:5][cH:6][cH:7]1",
        denticity=1,
        charge=-1,
        aliases=("anilido", "phenylamido"),
        description="Anilide",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    # Other common anionic
    "SCN": LigandInfo(
        smiles="N#C[S-]",
        mapped_smiles="[N:1]#[C:2][S-:3]",
        denticity=1,
        charge=-1,
        aliases=("thiocyanato", "thiocyanate"),
        description="Thiocyanate (S-bound)",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "NCS": LigandInfo(
        smiles="",
        mapped_smiles="",
        denticity=1,
        charge=-1,
        aliases=("isothiocyanato", "isothiocyanate"),
        description="Isothiocyanate (N-bound)",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "N3": LigandInfo(
        smiles="[N-]=[N+]=[N-]",
        mapped_smiles="[N-:1]=[N+:2]=[N-:3]",
        denticity=1,
        charge=-1,
        aliases=("azido", "azide"),
        description="Azide",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "NO2": LigandInfo(
        smiles="",
        mapped_smiles="",
        denticity=1,
        charge=-1,
        aliases=("nitrito", "nitrite"),
        description="Nitrite (N-bound nitro)",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "ONO": LigandInfo(
        smiles="O=N[O-]",
        mapped_smiles="[O:1]=[N:2][O-:3]",
        denticity=1,
        charge=-1,
        aliases=("nitrito-O",),
        description="Nitrite (O-bound nitrito)",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "SH": LigandInfo(
        smiles="[SH-]",
        mapped_smiles="[SH-:1]",
        denticity=1,
        charge=-1,
        aliases=("mercapto", "sulfhydryl", "thiolato"),
        description="Hydrosulfide/Thiolate",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "SPh": LigandInfo(
        smiles="[S-]c1ccccc1",
        mapped_smiles="[S-:1][c:2]1[cH:3][cH:4][cH:5][cH:6][cH:7]1",
        denticity=1,
        charge=-1,
        aliases=("thiophenolate", "phenylthiolato"),
        description="Thiophenolate",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "SMe": LigandInfo(
        smiles="C[S-]",
        mapped_smiles="[CH3:1][S-:2]",
        denticity=1,
        charge=-1,
        aliases=("methylthiolate", "methanethiolato"),
        description="Methylthiolate",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "StBu": LigandInfo(
        smiles="CC(C)(C)[S-]",
        mapped_smiles="[CH3:1][C:2]([CH3:3])([CH3:4])[S-:5]",
        denticity=1,
        charge=-1,
        aliases=("tert-butylthiolate",),
        description="tert-Butylthiolate",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "OCN": LigandInfo(
        smiles="N#C[O-]",
        mapped_smiles="[N:1]#[C:2][O-:3]",
        denticity=1,
        charge=-1,
        aliases=("cyanato", "cyanate"),
        description="Cyanate",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "NCO": LigandInfo(
        smiles="[N-]=C=O",
        mapped_smiles="[N-:1]=[C:2]=[O:3]",
        denticity=1,
        charge=-1,
        aliases=("isocyanato", "isocyanate"),
        description="Isocyanate",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    # Carboxylates
    "OBz": LigandInfo(
        smiles="O=C([O-])c1ccccc1",
        mapped_smiles="[O:1]=[C:2]([O-:3])[c:4]1[cH:5][cH:6][cH:7][cH:8][cH:9]1",
        denticity=1,
        charge=-1,
        aliases=("benzoato", "benzoate"),
        description="Benzoate",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "OPiv": LigandInfo(
        smiles="CC(C)(C)C(=O)[O-]",
        mapped_smiles="[CH3:1][C:2]([CH3:3])([CH3:4])[C:5](=[O:6])[O-:7]",
        denticity=1,
        charge=-1,
        aliases=("pivalato", "pivalate", "trimethylacetate"),
        description="Pivalate",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "O2CCF3": LigandInfo(
        smiles="O=C([O-])C(F)(F)F",
        mapped_smiles="[O:1]=[C:2]([O-:3])[C:4]([F:5])([F:6])[F:7]",
        denticity=1,
        charge=-1,
        aliases=("trifluoroacetato", "trifluoroacetate", "TFA"),
        description="Trifluoroacetate",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "formate": LigandInfo(
        smiles="O=C[O-]",
        mapped_smiles="[O:1]=[CH:2][O-:3]",
        denticity=1,
        charge=-1,
        aliases=("formato", "HCO2"),
        description="Formate",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    # Oxo and related
    "O": LigandInfo(
        smiles="[O-2]",
        mapped_smiles="[O-2:1]",
        denticity=1,
        charge=-2,
        aliases=("oxo", "oxide", "oxido"),
        description="Oxo",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "S": LigandInfo(
        smiles="[S-2]",
        mapped_smiles="[S-2:1]",
        denticity=1,
        charge=-2,
        aliases=("sulfido", "sulfide"),
        description="Sulfido",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "Se": LigandInfo(
        smiles="[Se-2]",
        mapped_smiles="[Se-2:1]",
        denticity=1,
        charge=-2,
        aliases=("selenido", "selenide"),
        description="Selenido",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "Te": LigandInfo(
        smiles="[Te-2]",
        mapped_smiles="[Te-2:1]",
        denticity=1,
        charge=-2,
        aliases=("tellurido", "telluride"),
        description="Tellurido",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "NR": LigandInfo(
        smiles="",
        mapped_smiles="",
        denticity=1,
        charge=-2,
        aliases=("imido",),
        description="Imido (generic)",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "NAr": LigandInfo(
        smiles="",
        mapped_smiles="",
        denticity=1,
        charge=-2,
        aliases=("arylimido",),
        description="Aryl imido",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "NtBu": LigandInfo(
        smiles="",
        mapped_smiles="",
        denticity=1,
        charge=-2,
        aliases=("tert-butylimido",),
        description="tert-Butyl imido",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "NAd": LigandInfo(
        smiles="",
        mapped_smiles="",
        denticity=1,
        charge=-2,
        aliases=("adamantylimido",),
        description="Adamantyl imido",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    # Borohydrides
    "BH4": LigandInfo(
        smiles="[BH4-]",
        mapped_smiles="[BH4-:1]",
        denticity=1,
        charge=-1,
        aliases=("borohydride", "tetrahydroborate"),
        description="Borohydride",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    # -------------------------------------------------------------------------
    # Bidentate Neutral Ligands
    # -------------------------------------------------------------------------
    "bpy": LigandInfo(
        smiles="c1ccc(-c2ccccn2)nc1",
        mapped_smiles="[cH:1]1[cH:2][cH:3][c:4](-[c:5]2[cH:6][cH:7][cH:8][cH:9][n:10]2)[n:11][cH:12]1",
        denticity=2,
        charge=0,
        aliases=("2,2'-bipyridine", "bipyridine", "bipy"),
        description="2,2'-Bipyridine",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "dtbbpy": LigandInfo(
        smiles="CC(C)(C)c1ccnc(-c2cc(C(C)(C)C)ccn2)c1",
        mapped_smiles="[CH3:1][C:2]([CH3:3])([CH3:4])[c:5]1[cH:6][cH:7][n:8][c:9](-[c:10]2[cH:11][c:12]([C:13]([CH3:14])([CH3:15])[CH3:16])[cH:17][cH:18][n:19]2)[cH:20]1",
        denticity=2,
        charge=0,
        aliases=(
            "4,4'-di-tert-butyl-2,2'-bipyridine",
            "di-tert-butylbipyridine",
            "dtbpy",
        ),
        description="4,4'-Di-tert-butyl-2,2'-bipyridine",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "phen": LigandInfo(
        smiles="c1cnc2c(c1)ccc1cccnc12",
        mapped_smiles="[cH:1]1[cH:2][n:3][c:4]2[c:5]([cH:6]1)[cH:7][cH:8][c:9]1[cH:10][cH:11][cH:12][n:13][c:14]21",
        denticity=2,
        charge=0,
        aliases=("1,10-phenanthroline", "phenanthroline"),
        description="1,10-Phenanthroline",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "en": LigandInfo(
        smiles="NCCN",
        mapped_smiles="[NH2:1][CH2:2][CH2:3][NH2:4]",
        denticity=2,
        charge=0,
        aliases=("ethylenediamine",),
        description="Ethylenediamine",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "cod": LigandInfo(
        smiles="C1=CCCC=CCC1",
        mapped_smiles="[CH:1]1=[CH:2][CH2:3][CH2:4][CH:5]=[CH:6][CH2:7][CH2:8]1",
        denticity=2,
        charge=0,
        aliases=("1,5-cyclooctadiene", "cyclooctadiene", "COD"),
        description="1,5-Cyclooctadiene (η⁴)",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "nbd": LigandInfo(
        smiles="C1=CC2C=CC1C2",
        mapped_smiles="[CH:1]1=[CH:2][CH:3]2[CH:4]=[CH:5][CH:6]1[CH2:7]2",
        denticity=2,
        charge=0,
        aliases=("norbornadiene", "2,5-norbornadiene"),
        description="Norbornadiene",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "dppe": LigandInfo(
        smiles="c1ccc(P(CCP(c2ccccc2)c2ccccc2)c2ccccc2)cc1",
        mapped_smiles="[cH:1]1[cH:2][cH:3][c:4]([P:5]([CH2:6][CH2:7][P:8]([c:9]2[cH:10][cH:11][cH:12][cH:13][cH:14]2)[c:15]2[cH:16][cH:17][cH:18][cH:19][cH:20]2)[c:21]2[cH:22][cH:23][cH:24][cH:25][cH:26]2)[cH:27][cH:28]1",
        denticity=2,
        charge=0,
        aliases=("1,2-bis(diphenylphosphino)ethane",),
        description="1,2-Bis(diphenylphosphino)ethane",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "dppm": LigandInfo(
        smiles="c1ccc(P(CP(c2ccccc2)c2ccccc2)c2ccccc2)cc1",
        mapped_smiles="[cH:1]1[cH:2][cH:3][c:4]([P:5]([CH2:6][P:7]([c:8]2[cH:9][cH:10][cH:11][cH:12][cH:13]2)[c:14]2[cH:15][cH:16][cH:17][cH:18][cH:19]2)[c:20]2[cH:21][cH:22][cH:23][cH:24][cH:25]2)[cH:26][cH:27]1",
        denticity=2,
        charge=0,
        aliases=("bis(diphenylphosphino)methane",),
        description="Bis(diphenylphosphino)methane",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    # Bipyridines and derivatives
    "4,4'-dmbpy": LigandInfo(
        smiles="Cc1ccnc(-c2cc(C)ccn2)c1",
        mapped_smiles="[CH3:1][c:2]1[cH:3][cH:4][n:5][c:6](-[c:7]2[cH:8][c:9]([CH3:10])[cH:11][cH:12][n:13]2)[cH:14]1",
        denticity=2,
        charge=0,
        aliases=("4,4'-dimethyl-2,2'-bipyridine", "dmb"),
        description="4,4'-Dimethyl-2,2'-bipyridine",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "5,5'-dmbpy": LigandInfo(
        smiles="Cc1ccc(-c2ccc(C)cn2)nc1",
        mapped_smiles="[CH3:1][c:2]1[cH:3][cH:4][c:5](-[c:6]2[cH:7][cH:8][c:9]([CH3:10])[cH:11][n:12]2)[n:13][cH:14]1",
        denticity=2,
        charge=0,
        aliases=("5,5'-dimethyl-2,2'-bipyridine",),
        description="5,5'-Dimethyl-2,2'-bipyridine",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "dCbpy": LigandInfo(
        smiles="O=C([O-])c1ccnc(-c2cc(C(=O)[O-])ccn2)c1",
        mapped_smiles="[O:1]=[C:2]([O-:3])[c:4]1[cH:5][cH:6][n:7][c:8](-[c:9]2[cH:10][c:11]([C:12](=[O:13])[O-:14])[cH:15][cH:16][n:17]2)[cH:18]1",
        denticity=2,
        charge=0,
        aliases=("4,4'-dicarboxy-2,2'-bipyridine",),
        description="4,4'-Dicarboxy-2,2'-bipyridine",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "dCEbpy": LigandInfo(
        smiles="CCOC(=O)c1ccnc(-c2cc(C(=O)OCC)ccn2)c1",
        mapped_smiles="[CH3:1][CH2:2][O:3][C:4](=[O:5])[c:6]1[cH:7][cH:8][n:9][c:10](-[c:11]2[cH:12][c:13]([C:14](=[O:15])[O:16][CH2:17][CH3:18])[cH:19][cH:20][n:21]2)[cH:22]1",
        denticity=2,
        charge=0,
        aliases=("4,4'-dicarboxyethyl-2,2'-bipyridine",),
        description="4,4'-Diethoxycarbonyl-2,2'-bipyridine",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    # Phenanthroline derivatives
    "dmp": LigandInfo(
        smiles="CC(=O)OI1(OC(C)=O)(OC(C)=O)OC(=O)c2ccccc21",
        mapped_smiles="[CH3:1][C:2](=[O:3])[O:4][I:5]1([O:6][C:7]([CH3:8])=[O:9])([O:10][C:11]([CH3:12])=[O:13])[O:14][C:15](=[O:16])[c:17]2[cH:18][cH:19][cH:20][cH:21][c:22]21",
        denticity=2,
        charge=0,
        aliases=("2,9-dimethyl-1,10-phenanthroline", "neocuproine"),
        description="2,9-Dimethyl-1,10-phenanthroline",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "dpp": LigandInfo(
        smiles="c1ccc(-c2ccc3ccc4ccc(-c5ccccc5)nc4c3n2)cc1",
        mapped_smiles="[cH:1]1[cH:2][cH:3][c:4](-[c:5]2[cH:6][cH:7][c:8]3[cH:9][cH:10][c:11]4[cH:12][cH:13][c:14](-[c:15]5[cH:16][cH:17][cH:18][cH:19][cH:20]5)[n:21][c:22]4[c:23]3[n:24]2)[cH:25][cH:26]1",
        denticity=2,
        charge=0,
        aliases=("2,9-diphenyl-1,10-phenanthroline", "bathocuproine"),
        description="2,9-Diphenyl-1,10-phenanthroline",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "tmp": LigandInfo(
        smiles="Cc1cnc2c(ccc3c(C)c(C)cnc32)c1C",
        mapped_smiles="[CH3:1][c:2]1[cH:3][n:4][c:5]2[c:6]([cH:7][cH:8][c:9]3[c:10]([CH3:11])[c:12]([CH3:13])[cH:14][n:15][c:16]23)[c:17]1[CH3:18]",
        denticity=2,
        charge=0,
        aliases=("3,4,7,8-tetramethyl-1,10-phenanthroline",),
        description="3,4,7,8-Tetramethyl-1,10-phenanthroline",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    # Bipyrimidines and related
    "bpm": LigandInfo(
        smiles="c1cnc(-c2ncccn2)nc1",
        mapped_smiles="[cH:1]1[cH:2][n:3][c:4](-[c:5]2[n:6][cH:7][cH:8][cH:9][n:10]2)[n:11][cH:12]1",
        denticity=2,
        charge=0,
        aliases=("2,2'-bipyrimidine", "bipyrimidine"),
        description="2,2'-Bipyrimidine",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "bpz": LigandInfo(
        smiles="c1cnc(-c2cnccn2)cn1",
        mapped_smiles="[cH:1]1[cH:2][n:3][c:4](-[c:5]2[cH:6][n:7][cH:8][cH:9][n:10]2)[cH:11][n:12]1",
        denticity=2,
        charge=0,
        aliases=("2,2'-bipyrazine", "bipyrazine"),
        description="2,2'-Bipyrazine",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    # Diphosphines (very important in catalysis)
    "dppp": LigandInfo(
        smiles="c1ccc(P(CCCP(c2ccccc2)c2ccccc2)c2ccccc2)cc1",
        mapped_smiles="[cH:1]1[cH:2][cH:3][c:4]([P:5]([CH2:6][CH2:7][CH2:8][P:9]([c:10]2[cH:11][cH:12][cH:13][cH:14][cH:15]2)[c:16]2[cH:17][cH:18][cH:19][cH:20][cH:21]2)[c:22]2[cH:23][cH:24][cH:25][cH:26][cH:27]2)[cH:28][cH:29]1",
        denticity=2,
        charge=0,
        aliases=("1,3-bis(diphenylphosphino)propane",),
        description="1,3-Bis(diphenylphosphino)propane",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "dppb": LigandInfo(
        smiles="c1ccc(P(CCCCP(c2ccccc2)c2ccccc2)c2ccccc2)cc1",
        mapped_smiles="[cH:1]1[cH:2][cH:3][c:4]([P:5]([CH2:6][CH2:7][CH2:8][CH2:9][P:10]([c:11]2[cH:12][cH:13][cH:14][cH:15][cH:16]2)[c:17]2[cH:18][cH:19][cH:20][cH:21][cH:22]2)[c:23]2[cH:24][cH:25][cH:26][cH:27][cH:28]2)[cH:29][cH:30]1",
        denticity=2,
        charge=0,
        aliases=("1,4-bis(diphenylphosphino)butane",),
        description="1,4-Bis(diphenylphosphino)butane",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "dppf": LigandInfo(
        smiles="c1ccc(P(c2ccccc2)[C]23[CH]4[CH]5[CH]6[CH]2[Fe]56432789[CH]3[CH]2[CH]7[C]8(P(c2ccccc2)c2ccccc2)[CH]39)cc1",
        mapped_smiles="[cH:1]1[cH:2][cH:3][c:4]([P:5]([c:6]2[cH:7][cH:8][cH:9][cH:10][cH:11]2)[C:12]23[CH:13]4[CH:14]5[CH:15]6[CH:16]2[Fe:17]34562789[CH:18]3[CH:19]2[CH:20]7[C:21]8([P:22]([c:23]2[cH:24][cH:25][cH:26][cH:27][cH:28]2)[c:29]2[cH:30][cH:31][cH:32][cH:33][cH:34]2)[CH:35]93)[cH:36][cH:37]1",
        denticity=2,
        charge=0,
        aliases=("1,1'-bis(diphenylphosphino)ferrocene",),
        description="1,1'-Bis(diphenylphosphino)ferrocene",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "dcpe": LigandInfo(
        smiles="C1CCC(P(CCP(C2CCCCC2)C2CCCCC2)C2CCCCC2)CC1",
        mapped_smiles="[CH2:1]1[CH2:2][CH2:3][CH:4]([P:5]([CH2:6][CH2:7][P:8]([CH:9]2[CH2:10][CH2:11][CH2:12][CH2:13][CH2:14]2)[CH:15]2[CH2:16][CH2:17][CH2:18][CH2:19][CH2:20]2)[CH:21]2[CH2:22][CH2:23][CH2:24][CH2:25][CH2:26]2)[CH2:27][CH2:28]1",
        denticity=2,
        charge=0,
        aliases=("1,2-bis(dicyclohexylphosphino)ethane",),
        description="1,2-Bis(dicyclohexylphosphino)ethane",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "dcpm": LigandInfo(
        smiles="C1CCC(P(CP(C2CCCCC2)C2CCCCC2)C2CCCCC2)CC1",
        mapped_smiles="[CH2:1]1[CH2:2][CH2:3][CH:4]([P:5]([CH2:6][P:7]([CH:8]2[CH2:9][CH2:10][CH2:11][CH2:12][CH2:13]2)[CH:14]2[CH2:15][CH2:16][CH2:17][CH2:18][CH2:19]2)[CH:20]2[CH2:21][CH2:22][CH2:23][CH2:24][CH2:25]2)[CH2:26][CH2:27]1",
        denticity=2,
        charge=0,
        aliases=("bis(dicyclohexylphosphino)methane",),
        description="Bis(dicyclohexylphosphino)methane",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "dmpe": LigandInfo(
        smiles="CP(C)CCP(C)C",
        mapped_smiles="[CH3:1][P:2]([CH3:3])[CH2:4][CH2:5][P:6]([CH3:7])[CH3:8]",
        denticity=2,
        charge=0,
        aliases=("1,2-bis(dimethylphosphino)ethane",),
        description="1,2-Bis(dimethylphosphino)ethane",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "depe": LigandInfo(
        smiles="CCP(CC)CCP(CC)CC",
        mapped_smiles="[CH3:1][CH2:2][P:3]([CH2:4][CH3:5])[CH2:6][CH2:7][P:8]([CH2:9][CH3:10])[CH2:11][CH3:12]",
        denticity=2,
        charge=0,
        aliases=("1,2-bis(diethylphosphino)ethane",),
        description="1,2-Bis(diethylphosphino)ethane",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "dippe": LigandInfo(
        smiles="CC(C)P(CCP(C(C)C)C(C)C)C(C)C",
        mapped_smiles="[CH3:1][CH:2]([CH3:3])[P:4]([CH2:5][CH2:6][P:7]([CH:8]([CH3:9])[CH3:10])[CH:11]([CH3:12])[CH3:13])[CH:14]([CH3:15])[CH3:16]",
        denticity=2,
        charge=0,
        aliases=("1,2-bis(diisopropylphosphino)ethane",),
        description="1,2-Bis(diisopropylphosphino)ethane",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "dtbpe": LigandInfo(
        smiles="CC(C)(C)P(CCP(C(C)(C)C)C(C)(C)C)C(C)(C)C",
        mapped_smiles="[CH3:1][C:2]([CH3:3])([CH3:4])[P:5]([CH2:6][CH2:7][P:8]([C:9]([CH3:10])([CH3:11])[CH3:12])[C:13]([CH3:14])([CH3:15])[CH3:16])[C:17]([CH3:18])([CH3:19])[CH3:20]",
        denticity=2,
        charge=0,
        aliases=("1,2-bis(di-tert-butylphosphino)ethane",),
        description="1,2-Bis(di-tert-butylphosphino)ethane",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    # BINAP and related chiral phosphines (critical in asymmetric catalysis)
    "BINAP": LigandInfo(
        smiles="c1ccc(P(c2ccccc2)c2ccc3ccccc3c2-c2c(P(c3ccccc3)c3ccccc3)ccc3ccccc23)cc1",
        mapped_smiles="[cH:1]1[cH:2][cH:3][c:4]([P:5]([c:6]2[cH:7][cH:8][cH:9][cH:10][cH:11]2)[c:12]2[cH:13][cH:14][c:15]3[cH:16][cH:17][cH:18][cH:19][c:20]3[c:21]2-[c:22]2[c:23]([P:24]([c:25]3[cH:26][cH:27][cH:28][cH:29][cH:30]3)[c:31]3[cH:32][cH:33][cH:34][cH:35][cH:36]3)[cH:37][cH:38][c:39]3[cH:40][cH:41][cH:42][cH:43][c:44]23)[cH:45][cH:46]1",
        denticity=2,
        charge=0,
        aliases=("2,2'-bis(diphenylphosphino)-1,1'-binaphthyl",),
        description="BINAP",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "TolBINAP": LigandInfo(
        smiles="Cc1ccc(P(c2ccc(C)cc2)c2ccc3ccccc3c2-c2c(P(c3ccc(C)cc3)c3ccc(C)cc3)ccc3ccccc23)cc1",
        mapped_smiles="[CH3:1][c:2]1[cH:3][cH:4][c:5]([P:6]([c:7]2[cH:8][cH:9][c:10]([CH3:11])[cH:12][cH:13]2)[c:14]2[cH:15][cH:16][c:17]3[cH:18][cH:19][cH:20][cH:21][c:22]3[c:23]2-[c:24]2[c:25]([P:26]([c:27]3[cH:28][cH:29][c:30]([CH3:31])[cH:32][cH:33]3)[c:34]3[cH:35][cH:36][c:37]([CH3:38])[cH:39][cH:40]3)[cH:41][cH:42][c:43]3[cH:44][cH:45][cH:46][cH:47][c:48]23)[cH:49][cH:50]1",
        denticity=2,
        charge=0,
        aliases=("2,2'-bis(di-p-tolylphosphino)-1,1'-binaphthyl",),
        description="Tol-BINAP",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "SEGPHOS": LigandInfo(
        smiles="c1ccc(P(c2ccccc2)c2ccc3c(c2-c2c(P(c4ccccc4)c4ccccc4)ccc4c2OCO4)OCO3)cc1",
        mapped_smiles="[cH:1]1[cH:2][cH:3][c:4]([P:5]([c:6]2[cH:7][cH:8][cH:9][cH:10][cH:11]2)[c:12]2[cH:13][cH:14][c:15]3[c:16]([c:17]2-[c:18]2[c:19]([P:20]([c:21]4[cH:22][cH:23][cH:24][cH:25][cH:26]4)[c:27]4[cH:28][cH:29][cH:30][cH:31][cH:32]4)[cH:33][cH:34][c:35]4[c:36]2[O:37][CH2:38][O:39]4)[O:40][CH2:41][O:42]3)[cH:43][cH:44]1",
        denticity=2,
        charge=0,
        aliases=("5,5'-bis(diphenylphosphino)-4,4'-bi-1,3-benzodioxole",),
        description="SEGPHOS",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "DM-SEGPHOS": LigandInfo(
        smiles="Cc1cc(C)cc(P(c2cc(C)cc(C)c2)c2ccc3c(c2-c2c(P(c4cc(C)cc(C)c4)c4cc(C)cc(C)c4)ccc4c2OCO4)OCO3)c1",
        mapped_smiles="[CH3:1][c:2]1[cH:3][c:4]([CH3:5])[cH:6][c:7]([P:8]([c:9]2[cH:10][c:11]([CH3:12])[cH:13][c:14]([CH3:15])[cH:16]2)[c:17]2[cH:18][cH:19][c:20]3[c:21]([c:22]2-[c:23]2[c:24]([P:25]([c:26]4[cH:27][c:28]([CH3:29])[cH:30][c:31]([CH3:32])[cH:33]4)[c:34]4[cH:35][c:36]([CH3:37])[cH:38][c:39]([CH3:40])[cH:41]4)[cH:42][cH:43][c:44]4[c:45]2[O:46][CH2:47][O:48]4)[O:49][CH2:50][O:51]3)[cH:52]1",
        denticity=2,
        charge=0,
        aliases=("5,5'-bis(di(3,5-xylyl)phosphino)-4,4'-bi-1,3-benzodioxole",),
        description="DM-SEGPHOS",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "DIFLUORPHOS": LigandInfo(
        smiles="FC1(F)Oc2ccc(P(c3ccccc3)c3ccccc3)c(-c3c(P(c4ccccc4)c4ccccc4)ccc4c3OC(F)(F)O4)c2O1",
        mapped_smiles="[F:1][C:2]1([F:3])[O:4][c:5]2[cH:6][cH:7][c:8]([P:9]([c:10]3[cH:11][cH:12][cH:13][cH:14][cH:15]3)[c:16]3[cH:17][cH:18][cH:19][cH:20][cH:21]3)[c:22](-[c:23]3[c:24]([P:25]([c:26]4[cH:27][cH:28][cH:29][cH:30][cH:31]4)[c:32]4[cH:33][cH:34][cH:35][cH:36][cH:37]4)[cH:38][cH:39][c:40]4[c:41]3[O:42][C:43]([F:44])([F:45])[O:46]4)[c:47]2[O:48]1",
        denticity=2,
        charge=0,
        aliases=(),
        description="DIFLUORPHOS",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "MeO-BIPHEP": LigandInfo(
        smiles="COc1cccc(P(c2ccccc2)c2ccccc2)c1-c1c(OC)cccc1P(c1ccccc1)c1ccccc1",
        mapped_smiles="[CH3:1][O:2][c:3]1[cH:4][cH:5][cH:6][c:7]([P:8]([c:9]2[cH:10][cH:11][cH:12][cH:13][cH:14]2)[c:15]2[cH:16][cH:17][cH:18][cH:19][cH:20]2)[c:21]1-[c:22]1[c:23]([O:24][CH3:25])[cH:26][cH:27][cH:28][c:29]1[P:30]([c:31]1[cH:32][cH:33][cH:34][cH:35][cH:36]1)[c:37]1[cH:38][cH:39][cH:40][cH:41][cH:42]1",
        denticity=2,
        charge=0,
        aliases=("6,6'-dimethoxy-2,2'-bis(diphenylphosphino)-1,1'-biphenyl",),
        description="MeO-BIPHEP",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    # Josiphos-type ligands
    "Josiphos": LigandInfo(
        smiles="C[C@@H]([C]1[CH][CH][CH][C]1P(c1ccccc1)c1ccccc1)P(C1CCCCC1)C1CCCCC1.[CH]1[CH][CH][CH][CH]1.[Fe]",
        mapped_smiles="[CH3:1][C@@H:2]([C:3]1[CH:4][CH:5][CH:6][C:7]1[P:8]([c:9]1[cH:10][cH:11][cH:12][cH:13][cH:14]1)[c:15]1[cH:16][cH:17][cH:18][cH:19][cH:20]1)[P:21]([CH:22]1[CH2:23][CH2:24][CH2:25][CH2:26][CH2:27]1)[CH:28]1[CH2:29][CH2:30][CH2:31][CH2:32][CH2:33]1.[CH:34]1[CH:35][CH:36][CH:37][CH:38]1.[Fe:39]",
        denticity=2,
        charge=0,
        aliases=(),
        description="Josiphos",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    # Other chiral ligands
    "CHIRAPHOS": LigandInfo(
        smiles="C[C@@H]([C@H](C)P(c1ccccc1)c1ccccc1)P(c1ccccc1)c1ccccc1",
        mapped_smiles="[CH3:1][C@@H:2]([C@H:3]([CH3:4])[P:5]([c:6]1[cH:7][cH:8][cH:9][cH:10][cH:11]1)[c:12]1[cH:13][cH:14][cH:15][cH:16][cH:17]1)[P:18]([c:19]1[cH:20][cH:21][cH:22][cH:23][cH:24]1)[c:25]1[cH:26][cH:27][cH:28][cH:29][cH:30]1",
        denticity=2,
        charge=0,
        aliases=("2,3-bis(diphenylphosphino)butane",),
        description="CHIRAPHOS",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "DIOP": LigandInfo(
        smiles="CC1(C)O[C@@H](CP(c2ccccc2)c2ccccc2)[C@H](CP(c2ccccc2)c2ccccc2)O1",
        mapped_smiles="[CH3:1][C:2]1([CH3:3])[O:4][C@@H:5]([CH2:6][P:7]([c:8]2[cH:9][cH:10][cH:11][cH:12][cH:13]2)[c:14]2[cH:15][cH:16][cH:17][cH:18][cH:19]2)[C@H:20]([CH2:21][P:22]([c:23]2[cH:24][cH:25][cH:26][cH:27][cH:28]2)[c:29]2[cH:30][cH:31][cH:32][cH:33][cH:34]2)[O:35]1",
        denticity=2,
        charge=0,
        aliases=(
            "2,3-O-isopropylidene-2,3-dihydroxy-1,4-bis(diphenylphosphino)butane",
        ),
        description="DIOP",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "DuPhos": LigandInfo(
        smiles="",
        mapped_smiles="",
        denticity=2,
        charge=0,
        aliases=("1,2-bis(phospholano)benzene",),
        description="DuPhos",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "BPE": LigandInfo(
        smiles="",
        mapped_smiles="",
        denticity=2,
        charge=0,
        aliases=("1,2-bis(phospholano)ethane",),
        description="BPE",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    # Diamines
    "tmeda": LigandInfo(
        smiles="CN(C)CCN(C)C",
        mapped_smiles="[CH3:1][N:2]([CH3:3])[CH2:4][CH2:5][N:6]([CH3:7])[CH3:8]",
        denticity=2,
        charge=0,
        aliases=("N,N,N',N'-tetramethylethylenediamine", "TMEDA"),
        description="TMEDA",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "dach": LigandInfo(
        smiles="NC1CCCCC1N",
        mapped_smiles="[NH2:1][CH:2]1[CH2:3][CH2:4][CH2:5][CH2:6][CH:7]1[NH2:8]",
        denticity=2,
        charge=0,
        aliases=("1,2-diaminocyclohexane", "chxn"),
        description="1,2-Diaminocyclohexane",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "dpen": LigandInfo(
        smiles="N[C@@H](c1ccccc1)[C@@H](N)c1ccccc1",
        mapped_smiles="[NH2:1][C@@H:2]([c:3]1[cH:4][cH:5][cH:6][cH:7][cH:8]1)[C@@H:9]([NH2:10])[c:11]1[cH:12][cH:13][cH:14][cH:15][cH:16]1",
        denticity=2,
        charge=0,
        aliases=("1,2-diphenylethylenediamine",),
        description="1,2-Diphenylethylenediamine",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "pn": LigandInfo(
        smiles="CC(N)CN",
        mapped_smiles="[CH3:1][CH:2]([NH2:3])[CH2:4][NH2:5]",
        denticity=2,
        charge=0,
        aliases=("1,2-diaminopropane", "propylenediamine"),
        description="1,2-Diaminopropane",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "bn": LigandInfo(
        smiles="CC(N)C(C)N",
        mapped_smiles="[CH3:1][CH:2]([NH2:3])[CH:4]([CH3:5])[NH2:6]",
        denticity=2,
        charge=0,
        aliases=("2,3-diaminobutane", "2,3-butanediamine"),
        description="2,3-Diaminobutane",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    # Diimine ligands
    "DAB": LigandInfo(
        smiles="N=CC=N",
        mapped_smiles="[NH:1]=[CH:2][CH:3]=[NH:4]",
        denticity=2,
        charge=0,
        aliases=("1,4-diazabutadiene",),
        description="1,4-Diazabutadiene",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "Ar-DAB": LigandInfo(
        smiles="",
        mapped_smiles="",
        denticity=2,
        charge=0,
        aliases=("N,N'-diaryl-1,4-diazabutadiene",),
        description="N,N'-Diaryl-1,4-diazabutadiene",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "dpp-BIAN": LigandInfo(
        smiles="CC(C)c1cccc(C(C)C)c1-c1cc2cccc3c2c(c1-c1c(C(C)C)cccc1C(C)C)C(=N)C3=N",
        mapped_smiles="[CH3:1][CH:2]([CH3:3])[c:4]1[cH:5][cH:6][cH:7][c:8]([CH:9]([CH3:10])[CH3:11])[c:12]1-[c:13]1[cH:14][c:15]2[cH:16][cH:17][cH:18][c:19]3[c:20]2[c:21]([c:22]1-[c:23]1[c:24]([CH:25]([CH3:26])[CH3:27])[cH:28][cH:29][cH:30][c:31]1[CH:32]([CH3:33])[CH3:34])[C:35](=[NH:36])[C:37]3=[NH:38]",
        denticity=2,
        charge=0,
        aliases=("bis(2,6-diisopropylphenyl)acenaphthenequinonediimine",),
        description="dpp-BIAN",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    # Mixed P,N donors
    "PHOX": LigandInfo(
        smiles="CC(C)(C)[C@H]1COC(c2ccccc2P(c2ccccc2)c2ccccc2)=N1",
        mapped_smiles="[CH3:1][C:2]([CH3:3])([CH3:4])[C@H:5]1[CH2:6][O:7][C:8]([c:9]2[cH:10][cH:11][cH:12][cH:13][c:14]2[P:15]([c:16]2[cH:17][cH:18][cH:19][cH:20][cH:21]2)[c:22]2[cH:23][cH:24][cH:25][cH:26][cH:27]2)=[N:28]1",
        denticity=2,
        charge=0,
        aliases=("phosphinooxazoline",),
        description="Phosphinooxazoline",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    # Schiff bases (common in coordination chemistry)
    "salen-H2": LigandInfo(
        smiles="Oc1ccccc1C=NCCN=Cc1ccccc1O",
        mapped_smiles="[OH:1][c:2]1[cH:3][cH:4][cH:5][cH:6][c:7]1[CH:8]=[N:9][CH2:10][CH2:11][N:12]=[CH:13][c:14]1[cH:15][cH:16][cH:17][cH:18][c:19]1[OH:20]",
        denticity=2,
        charge=0,
        aliases=("N,N'-bis(salicylidene)ethylenediamine",),
        description="Salen (neutral form)",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    # Other bidentate neutrals
    "dbm": LigandInfo(
        smiles="O=C(CC(=O)c1ccccc1)c1ccccc1",
        mapped_smiles="[O:1]=[C:2]([CH2:3][C:4](=[O:5])[c:6]1[cH:7][cH:8][cH:9][cH:10][cH:11]1)[c:12]1[cH:13][cH:14][cH:15][cH:16][cH:17]1",
        denticity=2,
        charge=0,
        aliases=("dibenzoylmethane",),
        description="Dibenzoylmethane (neutral)",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "dme": LigandInfo(
        smiles="COCCOC",
        mapped_smiles="[CH3:1][O:2][CH2:3][CH2:4][O:5][CH3:6]",
        denticity=2,
        charge=0,
        aliases=("1,2-dimethoxyethane", "glyme", "DME"),
        description="1,2-Dimethoxyethane",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "diglyme": LigandInfo(
        smiles="COCCOCCOC",
        mapped_smiles="[CH3:1][O:2][CH2:3][CH2:4][O:5][CH2:6][CH2:7][O:8][CH3:9]",
        denticity=2,
        charge=0,
        aliases=("diethylene glycol dimethyl ether",),
        description="Diglyme",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "OPPh3": LigandInfo(
        smiles="O=P(c1ccccc1)(c1ccccc1)c1ccccc1",
        mapped_smiles="[O:1]=[P:2]([c:3]1[cH:4][cH:5][cH:6][cH:7][cH:8]1)([c:9]1[cH:10][cH:11][cH:12][cH:13][cH:14]1)[c:15]1[cH:16][cH:17][cH:18][cH:19][cH:20]1",
        denticity=1,
        charge=0,
        aliases=("triphenylphosphine oxide",),
        description="Triphenylphosphine oxide",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    # -------------------------------------------------------------------------
    # Bidentate Anionic Ligands
    # -------------------------------------------------------------------------
    "acac": LigandInfo(
        smiles="CC(=O)[CH-]C(C)=O",
        mapped_smiles="[CH3:1][C:2](=[O:3])[CH-:4][C:5]([CH3:6])=[O:7]",
        denticity=2,
        charge=-1,
        aliases=("acetylacetonate", "acetylacetonato"),
        description="Acetylacetonate",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "ppy": LigandInfo(
        smiles="[c-]1ccccc1-c1ccccn1",
        mapped_smiles="[c-:1]1[cH:2][cH:3][cH:4][cH:5][c:6]1-[c:7]1[cH:8][cH:9][cH:10][cH:11][n:12]1",
        denticity=2,
        charge=-1,
        aliases=("2-phenylpyridine", "phenylpyridine", "phenylpyridinato"),
        description="2-Phenylpyridinate (C^N cyclometalating)",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "dfppy": LigandInfo(
        smiles="",
        mapped_smiles="",
        denticity=2,
        charge=-1,
        aliases=("2-(2,4-difluorophenyl)pyridine",),
        description="2-(2,4-Difluorophenyl)pyridinate",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "F2ppy": LigandInfo(
        smiles="",
        mapped_smiles="",
        denticity=2,
        charge=-1,
        aliases=("difluorophenylpyridine", "dFppy", "dF-ppy"),
        description="Difluorophenylpyridinate",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "dF(CF3)ppy": LigandInfo(
        smiles="",
        mapped_smiles="",
        denticity=2,
        charge=-1,
        aliases=("2-(2,4-Difluorophenyl)-5-(trifluoromethyl)pyridine", "dFCF3ppy"),
        description="2-(2,4-Difluorophenyl)-5-(trifluoromethyl)pyridine",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "pic": LigandInfo(
        smiles="O=C([O-])c1ccccn1",
        mapped_smiles="[O:1]=[C:2]([O-:3])[c:4]1[cH:5][cH:6][cH:7][cH:8][n:9]1",
        denticity=2,
        charge=-1,
        aliases=("picolinate", "picolinato"),
        description="Picolinate",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    # Cyclometalating C^N ligands
    "bzq": LigandInfo(
        smiles="c1ccc2c(c1)ccc1cccnc12",
        mapped_smiles="[cH:1]1[cH:2][cH:3][c:4]2[c:5]([cH:6]1)[cH:7][cH:8][c:9]1[cH:10][cH:11][cH:12][n:13][c:14]21",
        denticity=2,
        charge=-1,
        aliases=("benzo[h]quinoline", "7,8-benzoquinoline"),
        description="Benzo[h]quinolinate",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "thpy": LigandInfo(
        smiles="C1=C[S-]C(c2ccccn2)=C1",
        mapped_smiles="[CH:1]1=[CH:2][S-:3][C:4]([c:5]2[cH:6][cH:7][cH:8][cH:9][n:10]2)=[CH:11]1",
        denticity=2,
        charge=-1,
        aliases=("2-(2-thienyl)pyridine", "thienylpyridine"),
        description="2-(2-Thienyl)pyridinate",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "piq": LigandInfo(
        smiles="c1ccc(-c2nccc3ccccc23)cc1",
        mapped_smiles="[cH:1]1[cH:2][cH:3][c:4](-[c:5]2[n:6][cH:7][cH:8][c:9]3[cH:10][cH:11][cH:12][cH:13][c:14]23)[cH:15][cH:16]1",
        denticity=2,
        charge=-1,
        aliases=("1-phenylisoquinoline", "phenylisoquinoline"),
        description="1-Phenylisoquinolinate",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "mppy": LigandInfo(
        smiles="Cc1ccc(-c2ccccn2)cc1",
        mapped_smiles="[CH3:1][c:2]1[cH:3][cH:4][c:5](-[c:6]2[cH:7][cH:8][cH:9][cH:10][n:11]2)[cH:12][cH:13]1",
        denticity=2,
        charge=-1,
        aliases=("2-(4-methylphenyl)pyridine", "4-methyl-2-phenylpyridine"),
        description="2-(4-Methylphenyl)pyridinate",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "btp": LigandInfo(
        smiles="c1ccc(-c2cc3ccccc3s2)nc1",
        mapped_smiles="[cH:1]1[cH:2][cH:3][c:4](-[c:5]2[cH:6][c:7]3[cH:8][cH:9][cH:10][cH:11][c:12]3[s:13]2)[n:14][cH:15]1",
        denticity=2,
        charge=-1,
        aliases=("2-benzothienylpyridine",),
        description="2-Benzothienylpyridinate",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "pbpy": LigandInfo(
        smiles="c1ccc(-c2cccc(-c3ccccn3)n2)cc1",
        mapped_smiles="[cH:1]1[cH:2][cH:3][c:4](-[c:5]2[cH:6][cH:7][cH:8][c:9](-[c:10]3[cH:11][cH:12][cH:13][cH:14][n:15]3)[n:16]2)[cH:17][cH:18]1",
        denticity=2,
        charge=-1,
        aliases=("6-phenyl-2,2'-bipyridine",),
        description="6-Phenyl-2,2'-bipyridinate",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    # Beta-diketonates
    "hfac": LigandInfo(
        smiles="O=C(CC(=O)C(F)(F)F)C(F)(F)F",
        mapped_smiles="[O:1]=[C:2]([CH2:3][C:4](=[O:5])[C:6]([F:7])([F:8])[F:9])[C:10]([F:11])([F:12])[F:13]",
        denticity=2,
        charge=-1,
        aliases=("hexafluoroacetylacetonate", "1,1,1,5,5,5-hexafluoroacetylacetonate"),
        description="Hexafluoroacetylacetonate",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "tfac": LigandInfo(
        smiles="CC(=O)CC(=O)C(F)(F)F",
        mapped_smiles="[CH3:1][C:2](=[O:3])[CH2:4][C:5](=[O:6])[C:7]([F:8])([F:9])[F:10]",
        denticity=2,
        charge=-1,
        aliases=("trifluoroacetylacetonate",),
        description="Trifluoroacetylacetonate",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    # "dbm": LigandInfo(
    #     smiles="C(C1=CC=CC=C1)(=O)CC(C1=CC=CC=C1)=O",
    #     denticity=2,
    #     charge=-1,
    #     aliases=("dibenzoylmethanate",),
    #     description="Dibenzoylmethanate",
    # ),
    "thd": LigandInfo(
        smiles="CC(C)(C)C(=O)CC(=O)C(C)(C)C",
        mapped_smiles="[CH3:1][C:2]([CH3:3])([CH3:4])[C:5](=[O:6])[CH2:7][C:8](=[O:9])[C:10]([CH3:11])([CH3:12])[CH3:13]",
        denticity=2,
        charge=-1,
        aliases=("2,2,6,6-tetramethyl-3,5-heptanedionate", "tmhd", "dpm"),
        description="2,2,6,6-Tetramethyl-3,5-heptanedionate",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "fod": LigandInfo(
        smiles="CC(C)(C)C(=O)CC(=O)C(F)(F)C(F)(F)C(F)(F)F",
        mapped_smiles="[CH3:1][C:2]([CH3:3])([CH3:4])[C:5](=[O:6])[CH2:7][C:8](=[O:9])[C:10]([F:11])([F:12])[C:13]([F:14])([F:15])[C:16]([F:17])([F:18])[F:19]",
        denticity=2,
        charge=-1,
        aliases=("6,6,7,7,8,8,8-heptafluoro-2,2-dimethyl-3,5-octanedionate",),
        description="FOD",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "trop": LigandInfo(
        smiles="O=c1cccccc1O",
        mapped_smiles="[O:1]=[c:2]1[cH:3][cH:4][cH:5][cH:6][cH:7][c:8]1[OH:9]",
        denticity=2,
        charge=-1,
        aliases=("tropolonate",),
        description="Tropolonate",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    # Carboxylates (bridging/chelating)
    "OAc-bi": LigandInfo(
        smiles="",
        mapped_smiles="",
        denticity=2,
        charge=-1,
        aliases=("acetate-O,O'",),
        description="Acetate (bidentate)",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "CO3": LigandInfo(
        smiles="",
        mapped_smiles="",
        denticity=2,
        charge=-2,
        aliases=("carbonato", "carbonate"),
        description="Carbonate",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "SO4": LigandInfo(
        smiles="",
        mapped_smiles="",
        denticity=2,
        charge=-2,
        aliases=("sulfato", "sulfate"),
        description="Sulfate",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    # N,N chelates (anionic)
    "pz": LigandInfo(
        smiles="",
        mapped_smiles="",
        denticity=1,
        charge=-1,
        aliases=("pyrazolato", "pyrazolate"),
        description="Pyrazolate",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "pypz": LigandInfo(
        smiles="",
        mapped_smiles="",
        denticity=2,
        charge=-1,
        aliases=("3-(2-pyridyl)pyrazolate",),
        description="3-(2-Pyridyl)pyrazolate",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "indazolato": LigandInfo(
        smiles="",
        mapped_smiles="",
        denticity=2,
        charge=-1,
        aliases=("indazolate",),
        description="Indazolate",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    # N,O chelates
    "quinolinolate": LigandInfo(
        smiles="",
        mapped_smiles="",
        denticity=2,
        charge=-1,
        aliases=("8-hydroxyquinolinate", "oxinate", "Q"),
        description="8-Quinolinolate",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "glycinato": LigandInfo(
        smiles="",
        mapped_smiles="",
        denticity=2,
        charge=-1,
        aliases=("glycinate", "gly"),
        description="Glycinate",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "alaninato": LigandInfo(
        smiles="",
        mapped_smiles="",
        denticity=2,
        charge=-1,
        aliases=("alaninate", "ala"),
        description="Alaninate",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "salicylate": LigandInfo(
        smiles="",
        mapped_smiles="",
        denticity=2,
        charge=-1,
        aliases=("salicylato", "sal"),
        description="Salicylate",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "oxalato-mono": LigandInfo(
        smiles="",
        mapped_smiles="",
        denticity=2,
        charge=-1,
        aliases=(),
        description="Oxalate (monoanionic, monodentate)",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    # O,O chelates
    "catecholato": LigandInfo(
        smiles="",
        mapped_smiles="",
        denticity=2,
        charge=-2,
        aliases=("catecholate", "cat"),
        description="Catecholate",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "semiquinone": LigandInfo(
        smiles="",
        mapped_smiles="",
        denticity=2,
        charge=-1,
        aliases=("sq",),
        description="Semiquinone",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    # S,S chelates
    "S2CNMe2": LigandInfo(
        smiles="CN(C)C(=S)[S-]",
        mapped_smiles="[CH3:1][N:2]([CH3:3])[C:4](=[S:5])[S-:6]",
        denticity=2,
        charge=-1,
        aliases=("dimethyldithiocarbamate", "Me2dtc"),
        description="Dimethyldithiocarbamate",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "S2CNEt2": LigandInfo(
        smiles="CCN(CC)C(=S)[S-]",
        mapped_smiles="[CH3:1][CH2:2][N:3]([CH2:4][CH3:5])[C:6](=[S:7])[S-:8]",
        denticity=2,
        charge=-1,
        aliases=("diethyldithiocarbamate", "Et2dtc", "dtc"),
        description="Diethyldithiocarbamate",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "S2COEt": LigandInfo(
        smiles="CCOC(=S)[S-]",
        mapped_smiles="[CH3:1][CH2:2][O:3][C:4](=[S:5])[S-:6]",
        denticity=2,
        charge=-1,
        aliases=("ethylxanthate", "xanthate"),
        description="Ethyl xanthate",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "S2PPh2": LigandInfo(
        smiles="S=P([S-])(c1ccccc1)c1ccccc1",
        mapped_smiles="[S:1]=[P:2]([S-:3])([c:4]1[cH:5][cH:6][cH:7][cH:8][cH:9]1)[c:10]1[cH:11][cH:12][cH:13][cH:14][cH:15]1",
        denticity=2,
        charge=-1,
        aliases=("diphenylphosphinodithioate", "dtp"),
        description="Diphenylphosphinodithioate",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "bdt": LigandInfo(
        smiles="[S-]c1ccccc1[S-]",
        mapped_smiles="[S-:1][c:2]1[cH:3][cH:4][cH:5][cH:6][c:7]1[S-:8]",
        denticity=2,
        charge=-2,
        aliases=("1,2-benzenedithiolate", "benzene-1,2-dithiolate"),
        description="1,2-Benzenedithiolate",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "mnt": LigandInfo(
        smiles="",
        mapped_smiles="",
        denticity=2,
        charge=-2,
        aliases=("maleonitriledithiolate", "1,2-dicyanoethylene-1,2-dithiolate"),
        description="Maleonitriledithiolate",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "dmit": LigandInfo(
        smiles="S=c1sc([S-])c([S-])s1",
        mapped_smiles="[S:1]=[c:2]1[s:3][c:4]([S-:5])[c:6]([S-:7])[s:8]1",
        denticity=2,
        charge=-2,
        aliases=("2-thioxo-1,3-dithiole-4,5-dithiolate",),
        description="dmit",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "edt": LigandInfo(
        smiles="[S-]CC[S-]",
        mapped_smiles="[S-:1][CH2:2][CH2:3][S-:4]",
        denticity=2,
        charge=-2,
        aliases=("ethane-1,2-dithiolate", "1,2-ethanedithiolate"),
        description="Ethane-1,2-dithiolate",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "tdt": LigandInfo(
        smiles="Cc1ccc([S-])c([S-])c1",
        mapped_smiles="[CH3:1][c:2]1[cH:3][cH:4][c:5]([S-:6])[c:7]([S-:8])[cH:9]1",
        denticity=2,
        charge=-2,
        aliases=("toluene-3,4-dithiolate",),
        description="Toluene-3,4-dithiolate",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    # -------------------------------------------------------------------------
    # Tridentate Neutral Ligands
    # -------------------------------------------------------------------------
    # Terpyridines
    "tpy": LigandInfo(
        smiles="c1ccc(-c2cccc(-c3ccccn3)n2)nc1",
        mapped_smiles="[cH:1]1[cH:2][cH:3][c:4](-[c:5]2[cH:6][cH:7][cH:8][c:9](-[c:10]3[cH:11][cH:12][cH:13][cH:14][n:15]3)[n:16]2)[n:17][cH:18]1",
        denticity=3,
        charge=0,
        aliases=("2,2':6',2''-terpyridine", "terpyridine", "terpy"),
        description="2,2':6',2''-Terpyridine",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "ttpy": LigandInfo(
        smiles="Cc1ccc(-c2cc(-c3ccccn3)nc(-c3ccccn3)c2)cc1",
        mapped_smiles="[CH3:1][c:2]1[cH:3][cH:4][c:5](-[c:6]2[cH:7][c:8](-[c:9]3[cH:10][cH:11][cH:12][cH:13][n:14]3)[n:15][c:16](-[c:17]3[cH:18][cH:19][cH:20][cH:21][n:22]3)[cH:23]2)[cH:24][cH:25]1",
        denticity=3,
        charge=0,
        aliases=("4'-p-tolyl-2,2':6',2''-terpyridine",),
        description="4'-p-Tolyl-terpyridine",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "tBu3tpy": LigandInfo(
        smiles="CC(C)(C)c1ccnc(-c2cc(C(C)(C)C)cc(-c3cc(C(C)(C)C)ccn3)n2)c1",
        mapped_smiles="[CH3:1][C:2]([CH3:3])([CH3:4])[c:5]1[cH:6][cH:7][n:8][c:9](-[c:10]2[cH:11][c:12]([C:13]([CH3:14])([CH3:15])[CH3:16])[cH:17][c:18](-[c:19]3[cH:20][c:21]([C:22]([CH3:23])([CH3:24])[CH3:25])[cH:26][cH:27][n:28]3)[n:29]2)[cH:30]1",
        denticity=3,
        charge=0,
        aliases=("4,4',4''-tri-tert-butyl-2,2':6',2''-terpyridine",),
        description="4,4',4''-Tri-tert-butylterpyridine",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    # Pincer ligands (hugely important class)
    "PNP": LigandInfo(
        smiles="c1ccc(P(c2ccccc2)c2cccc(P(c3ccccc3)c3ccccc3)n2)cc1",
        mapped_smiles="[cH:1]1[cH:2][cH:3][c:4]([P:5]([c:6]2[cH:7][cH:8][cH:9][cH:10][cH:11]2)[c:12]2[cH:13][cH:14][cH:15][c:16]([P:17]([c:18]3[cH:19][cH:20][cH:21][cH:22][cH:23]3)[c:24]3[cH:25][cH:26][cH:27][cH:28][cH:29]3)[n:30]2)[cH:31][cH:32]1",
        denticity=3,
        charge=0,
        aliases=("bis(phosphino)pyridine",),
        description="PNP pincer",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "PCP": LigandInfo(
        smiles="",
        mapped_smiles="",
        denticity=3,
        charge=0,
        aliases=("bis(phosphino)aryl",),
        description="PCP pincer",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "NCN": LigandInfo(
        smiles="",
        mapped_smiles="",
        denticity=3,
        charge=0,
        aliases=("bis(amino)aryl",),
        description="NCN pincer",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "SCS": LigandInfo(
        smiles="",
        mapped_smiles="",
        denticity=3,
        charge=0,
        aliases=("bis(thio)aryl",),
        description="SCS pincer",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "CNC": LigandInfo(
        smiles="",
        mapped_smiles="",
        denticity=3,
        charge=0,
        aliases=("bis(NHC)pyridine",),
        description="CNC pincer (bis-NHC)",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    # PyBOX
    "PyBOX": LigandInfo(
        smiles="c1cc(C2=NCCO2)nc(C2=NCCO2)c1",
        mapped_smiles="[cH:1]1[cH:2][c:3]([C:4]2=[N:5][CH2:6][CH2:7][O:8]2)[n:9][c:10]([C:11]2=[N:12][CH2:13][CH2:14][O:15]2)[cH:16]1",
        denticity=3,
        charge=0,
        aliases=("pyridinebisoxazoline", "pybox"),
        description="2,6-Bis(oxazolinyl)pyridine",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "iPr-PyBOX": LigandInfo(
        smiles="CC(C)C1COC(c2cccc(C3=NC(C(C)C)CO3)n2)=N1",
        mapped_smiles="[CH3:1][CH:2]([CH3:3])[CH:4]1[CH2:5][O:6][C:7]([c:8]2[cH:9][cH:10][cH:11][c:12]([C:13]3=[N:14][CH:15]([CH:16]([CH3:17])[CH3:18])[CH2:19][O:20]3)[n:21]2)=[N:22]1",
        denticity=3,
        charge=0,
        aliases=(),
        description="2,6-Bis(4-isopropyl-2-oxazolinyl)pyridine",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "Ph-PyBOX": LigandInfo(
        smiles="c1ccc(C2COC(c3cccc(C4=NC(c5ccccc5)CO4)n3)=N2)cc1",
        mapped_smiles="[cH:1]1[cH:2][cH:3][c:4]([CH:5]2[CH2:6][O:7][C:8]([c:9]3[cH:10][cH:11][cH:12][c:13]([C:14]4=[N:15][CH:16]([c:17]5[cH:18][cH:19][cH:20][cH:21][cH:22]5)[CH2:23][O:24]4)[n:25]3)=[N:26]2)[cH:27][cH:28]1",
        denticity=3,
        charge=0,
        aliases=(),
        description="2,6-Bis(4-phenyl-2-oxazolinyl)pyridine",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    # Triphosphines
    "triphos": LigandInfo(
        smiles="CC(CP(c1ccccc1)c1ccccc1)(CP(c1ccccc1)c1ccccc1)CP(c1ccccc1)c1ccccc1",
        mapped_smiles="[CH3:1][C:2]([CH2:3][P:4]([c:5]1[cH:6][cH:7][cH:8][cH:9][cH:10]1)[c:11]1[cH:12][cH:13][cH:14][cH:15][cH:16]1)([CH2:17][P:18]([c:19]1[cH:20][cH:21][cH:22][cH:23][cH:24]1)[c:25]1[cH:26][cH:27][cH:28][cH:29][cH:30]1)[CH2:31][P:32]([c:33]1[cH:34][cH:35][cH:36][cH:37][cH:38]1)[c:39]1[cH:40][cH:41][cH:42][cH:43][cH:44]1",
        denticity=3,
        charge=0,
        aliases=("1,1,1-tris(diphenylphosphinomethyl)ethane", "MeC(CH2PPh2)3"),
        description="Triphos",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    # Triamines
    "dien": LigandInfo(
        smiles="NCCNCCN",
        mapped_smiles="[NH2:1][CH2:2][CH2:3][NH:4][CH2:5][CH2:6][NH2:7]",
        denticity=3,
        charge=0,
        aliases=("diethylenetriamine",),
        description="Diethylenetriamine",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "tacn": LigandInfo(
        smiles="C1CNCCNCCN1",
        mapped_smiles="[CH2:1]1[CH2:2][NH:3][CH2:4][CH2:5][NH:6][CH2:7][CH2:8][NH:9]1",
        denticity=3,
        charge=0,
        aliases=("1,4,7-triazacyclononane",),
        description="1,4,7-Triazacyclononane",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "Me3tacn": LigandInfo(
        smiles="CN1CCN(C)CCN(C)CC1",
        mapped_smiles="[CH3:1][N:2]1[CH2:3][CH2:4][N:5]([CH3:6])[CH2:7][CH2:8][N:9]([CH3:10])[CH2:11][CH2:12]1",
        denticity=3,
        charge=0,
        aliases=("1,4,7-trimethyl-1,4,7-triazacyclononane",),
        description="1,4,7-Trimethyl-1,4,7-triazacyclononane",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    # Scorpionate-type (neutral)
    "Tp": LigandInfo(
        smiles="",
        mapped_smiles="",
        denticity=3,
        charge=-1,
        aliases=("hydrotris(pyrazolyl)borate", "trispyrazolylborate"),
        description="Hydrotris(pyrazolyl)borate",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "Tp*": LigandInfo(
        smiles="",
        mapped_smiles="",
        denticity=3,
        charge=-1,
        aliases=("hydrotris(3,5-dimethylpyrazolyl)borate",),
        description="Hydrotris(3,5-dimethylpyrazolyl)borate",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "TpiPr2": LigandInfo(
        smiles="",
        mapped_smiles="",
        denticity=3,
        charge=-1,
        aliases=("hydrotris(3,5-diisopropylpyrazolyl)borate",),
        description="Hydrotris(3,5-diisopropylpyrazolyl)borate",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    # Other
    "bpa": LigandInfo(
        smiles="c1ccc(CNCc2ccccn2)nc1",
        mapped_smiles="[cH:1]1[cH:2][cH:3][c:4]([CH2:5][NH:6][CH2:7][c:8]2[cH:9][cH:10][cH:11][cH:12][n:13]2)[n:14][cH:15]1",
        denticity=3,
        charge=0,
        aliases=("bis(2-pyridylmethyl)amine",),
        description="Bis(2-pyridylmethyl)amine",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "bpea": LigandInfo(
        smiles="CCN(Cc1ccccn1)Cc1ccccn1",
        mapped_smiles="[CH3:1][CH2:2][N:3]([CH2:4][c:5]1[cH:6][cH:7][cH:8][cH:9][n:10]1)[CH2:11][c:12]1[cH:13][cH:14][cH:15][cH:16][n:17]1",
        denticity=3,
        charge=0,
        aliases=("N,N-bis(2-pyridylmethyl)ethylamine",),
        description="N,N-Bis(2-pyridylmethyl)ethylamine",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "dap": LigandInfo(
        smiles="CC(=O)c1cccc(C(C)=O)n1",
        mapped_smiles="[CH3:1][C:2](=[O:3])[c:4]1[cH:5][cH:6][cH:7][c:8]([C:9]([CH3:10])=[O:11])[n:12]1",
        denticity=3,
        charge=0,
        aliases=("2,6-Diacetylpyridine",),
        description="2,6-Diacetylpyridine",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    # -------------------------------------------------------------------------
    # Tridentate Anionic Ligands
    # -------------------------------------------------------------------------
    # Pincer anionic
    "PCP-": LigandInfo(
        smiles="",
        mapped_smiles="",
        denticity=3,
        charge=-1,
        aliases=(),
        description="PCP pincer (cyclometalated)",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "NCN-": LigandInfo(
        smiles="",
        mapped_smiles="",
        denticity=3,
        charge=-1,
        aliases=(),
        description="NCN pincer (cyclometalated)",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "PNP-": LigandInfo(
        smiles="",
        mapped_smiles="",
        denticity=3,
        charge=-1,
        aliases=(),
        description="PNP pincer (amido form)",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    # Bis(imino)pyridine
    "PDI": LigandInfo(
        smiles="",
        mapped_smiles="",
        denticity=3,
        charge=-1,
        aliases=("bis(imino)pyridine", "pyridinediimine"),
        description="Pyridine-2,6-diimine (reduced form)",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    # Corroles
    "corrole": LigandInfo(
        smiles="",
        mapped_smiles="",
        denticity=4,
        charge=-3,
        aliases=(),
        description="Corrole (tridentate in some counting)",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    # -------------------------------------------------------------------------
    # Tetradentate Neutral Ligands
    # -------------------------------------------------------------------------
    # Macrocycles
    "cyclam": LigandInfo(
        smiles="C1CNCCNCCCNCCNC1",
        mapped_smiles="[CH2:1]1[CH2:2][NH:3][CH2:4][CH2:5][NH:6][CH2:7][CH2:8][CH2:9][NH:10][CH2:11][CH2:12][NH:13][CH2:14]1",
        denticity=4,
        charge=0,
        aliases=("1,4,8,11-tetraazacyclotetradecane",),
        description="Cyclam",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "cyclen": LigandInfo(
        smiles="C1CNCCNCCNCCN1",
        mapped_smiles="[CH2:1]1[CH2:2][NH:3][CH2:4][CH2:5][NH:6][CH2:7][CH2:8][NH:9][CH2:10][CH2:11][NH:12]1",
        denticity=4,
        charge=0,
        aliases=("1,4,7,10-tetraazacyclododecane",),
        description="Cyclen",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "Me4cyclam": LigandInfo(
        smiles="CN1CCCN(C)CCN(C)CCCN(C)CC1",
        mapped_smiles="[CH3:1][N:2]1[CH2:3][CH2:4][CH2:5][N:6]([CH3:7])[CH2:8][CH2:9][N:10]([CH3:11])[CH2:12][CH2:13][CH2:14][N:15]([CH3:16])[CH2:17][CH2:18]1",
        denticity=4,
        charge=0,
        aliases=("1,4,8,11-tetramethyl-1,4,8,11-tetraazacyclotetradecane",),
        description="Tetramethylcyclam",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    # Linear tetradentate
    "trien": LigandInfo(
        smiles="NCCNCCNCCN",
        mapped_smiles="[NH2:1][CH2:2][CH2:3][NH:4][CH2:5][CH2:6][NH:7][CH2:8][CH2:9][NH2:10]",
        denticity=4,
        charge=0,
        aliases=("triethylenetetramine",),
        description="Triethylenetetramine",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    # Tetraphosphines
    "PP3": LigandInfo(
        smiles="",
        mapped_smiles="",
        denticity=4,
        charge=0,
        aliases=("tris(2-(diphenylphosphino)ethyl)phosphine",),
        description="PP3",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    # -------------------------------------------------------------------------
    # Tetradentate Anionic Ligands
    # -------------------------------------------------------------------------
    # Porphyrins (crucial ligand class)
    "TPP": LigandInfo(
        smiles="C1=Cc2nc1c(-c1ccccc1)c1ccc([nH]1)c(-c1ccccc1)c1nc(c(-c3ccccc3)c3ccc([nH]3)c2-c2ccccc2)C=C1",
        mapped_smiles="[CH:1]1=[CH:2][c:3]2[n:4][c:5]1[c:6](-[c:7]1[cH:8][cH:9][cH:10][cH:11][cH:12]1)[c:13]1[cH:14][cH:15][c:16]([nH:17]1)[c:18](-[c:19]1[cH:20][cH:21][cH:22][cH:23][cH:24]1)[c:25]1[n:26][c:27]([c:28](-[c:29]3[cH:30][cH:31][cH:32][cH:33][cH:34]3)[c:35]3[cH:36][cH:37][c:38]([nH:39]3)[c:40]2-[c:41]2[cH:42][cH:43][cH:44][cH:45][cH:46]2)[CH:47]=[CH:48]1",
        denticity=4,
        charge=-2,
        aliases=("tetraphenylporphyrin", "5,10,15,20-tetraphenylporphyrin"),
        description="Tetraphenylporphyrin",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "OEP": LigandInfo(
        smiles="CCC1=C(CC)c2cc3[nH]c(cc4nc(cc5[nH]c(cc1n2)c(CC)c5CC)C(CC)=C4CC)c(CC)c3CC",
        mapped_smiles="[CH3:1][CH2:2][C:3]1=[C:4]([CH2:5][CH3:6])[c:7]2[cH:8][c:9]3[nH:10][c:11]([cH:12][c:13]4[n:14][c:15]([cH:16][c:17]5[nH:18][c:19]([cH:20][c:21]1[n:22]2)[c:23]([CH2:24][CH3:25])[c:26]5[CH2:27][CH3:28])[C:29]([CH2:30][CH3:31])=[C:32]4[CH2:33][CH3:34])[c:35]([CH2:36][CH3:37])[c:38]3[CH2:39][CH3:40]",
        denticity=4,
        charge=-2,
        aliases=("octaethylporphyrin", "2,3,7,8,12,13,17,18-octaethylporphyrin"),
        description="Octaethylporphyrin",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "TMP": LigandInfo(
        smiles="Cc1cc(C)c(-c2c3nc(c(-c4c(C)cc(C)cc4C)c4ccc([nH]4)c(-c4c(C)cc(C)cc4C)c4nc(c(-c5c(C)cc(C)cc5C)c5ccc2[nH]5)C=C4)C=C3)c(C)c1",
        mapped_smiles="[CH3:1][c:2]1[cH:3][c:4]([CH3:5])[c:6](-[c:7]2[c:8]3[n:9][c:10]([c:11](-[c:12]4[c:13]([CH3:14])[cH:15][c:16]([CH3:17])[cH:18][c:19]4[CH3:20])[c:21]4[cH:22][cH:23][c:24]([nH:25]4)[c:26](-[c:27]4[c:28]([CH3:29])[cH:30][c:31]([CH3:32])[cH:33][c:34]4[CH3:35])[c:36]4[n:37][c:38]([c:39](-[c:40]5[c:41]([CH3:42])[cH:43][c:44]([CH3:45])[cH:46][c:47]5[CH3:48])[c:49]5[cH:50][cH:51][c:52]2[nH:53]5)[CH:54]=[CH:55]4)[CH:56]=[CH:57]3)[c:58]([CH3:59])[cH:60]1",
        denticity=4,
        charge=-2,
        aliases=("tetramesitylporphyrin",),
        description="Tetramesitylporphyrin",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "por": LigandInfo(
        smiles="C1=Cc2cc3ccc(cc4nc(cc5ccc(cc1n2)[nH]5)C=C4)[nH]3",
        mapped_smiles="[CH:1]1=[CH:2][c:3]2[cH:4][c:5]3[cH:6][cH:7][c:8]([cH:9][c:10]4[n:11][c:12]([cH:13][c:14]5[cH:15][cH:16][c:17]([cH:18][c:19]1[n:20]2)[nH:21]5)[CH:22]=[CH:23]4)[nH:24]3",
        denticity=4,
        charge=-2,
        aliases=("porphyrin", "porphyrinato"),
        description="Porphyrin (generic)",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "TPFPP": LigandInfo(
        smiles="",
        mapped_smiles="",
        denticity=4,
        charge=-2,
        aliases=("tetrakis(pentafluorophenyl)porphyrin",),
        description="Tetrakis(pentafluorophenyl)porphyrin",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    # Phthalocyanines
    "Pc": LigandInfo(
        smiles="",
        mapped_smiles="",
        denticity=4,
        charge=-2,
        aliases=("phthalocyanine", "phthalocyaninato"),
        description="Phthalocyanine",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    # Salen-type
    "salen": LigandInfo(
        smiles="",
        mapped_smiles="",
        denticity=4,
        charge=-2,
        aliases=("N,N'-ethylenebis(salicylideneiminato)",),
        description="Salen",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "salphen": LigandInfo(
        smiles="",
        mapped_smiles="",
        denticity=4,
        charge=-2,
        aliases=("N,N'-phenylenebis(salicylideneiminato)",),
        description="Salphen",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "salophen": LigandInfo(
        smiles="",
        mapped_smiles="",
        denticity=4,
        charge=-2,
        aliases=(),
        description="Salophen",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "salcn": LigandInfo(
        smiles="",
        mapped_smiles="",
        denticity=4,
        charge=-2,
        aliases=("N,N'-cyclohexanebis(salicylideneiminato)",),
        description="Salen-cyclohexanediamine",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "Jacobsen": LigandInfo(
        smiles="",
        mapped_smiles="",
        denticity=4,
        charge=-2,
        aliases=("Jacobsen's salen",),
        description="Jacobsen's salen ligand",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    # -------------------------------------------------------------------------
    # Penta/Hexadentate Neutral Ligands
    # -------------------------------------------------------------------------
    # Pentadentate
    "Me6tren": LigandInfo(
        smiles="CN(C)CCN(CCN(C)C)CCN(C)C",
        mapped_smiles="[CH3:1][N:2]([CH3:3])[CH2:4][CH2:5][N:6]([CH2:7][CH2:8][N:9]([CH3:10])[CH3:11])[CH2:12][CH2:13][N:14]([CH3:15])[CH3:16]",
        denticity=4,
        charge=0,
        aliases=("tris(2-dimethylaminoethyl)amine",),
        description="Tris(2-dimethylaminoethyl)amine",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "tpa": LigandInfo(
        smiles="c1ccc(CN(Cc2ccccn2)Cc2ccccn2)nc1",
        mapped_smiles="[cH:1]1[cH:2][cH:3][c:4]([CH2:5][N:6]([CH2:7][c:8]2[cH:9][cH:10][cH:11][cH:12][n:13]2)[CH2:14][c:15]2[cH:16][cH:17][cH:18][cH:19][n:20]2)[n:21][cH:22]1",
        denticity=4,
        charge=0,
        aliases=("tris(2-pyridylmethyl)amine",),
        description="Tris(2-pyridylmethyl)amine",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    # Hexadentate
    "EDTA": LigandInfo(
        smiles="O=C(O)CN(CCN(CC(=O)O)CC(=O)O)CC(=O)O",
        mapped_smiles="[O:1]=[C:2]([OH:3])[CH2:4][N:5]([CH2:6][CH2:7][N:8]([CH2:9][C:10](=[O:11])[OH:12])[CH2:13][C:14](=[O:15])[OH:16])[CH2:17][C:18](=[O:19])[OH:20]",
        denticity=6,
        charge=-4,
        aliases=("ethylenediaminetetraacetate", "edta"),
        description="Ethylenediaminetetraacetate",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "DTPA": LigandInfo(
        smiles="O=C(O)CN(CCN(CC(=O)O)CC(=O)O)CCN(CC(=O)O)CC(=O)O",
        mapped_smiles="[O:1]=[C:2]([OH:3])[CH2:4][N:5]([CH2:6][CH2:7][N:8]([CH2:9][C:10](=[O:11])[OH:12])[CH2:13][C:14](=[O:15])[OH:16])[CH2:17][CH2:18][N:19]([CH2:20][C:21](=[O:22])[OH:23])[CH2:24][C:25](=[O:26])[OH:27]",
        denticity=8,
        charge=-5,
        aliases=("diethylenetriaminepentaacetate",),
        description="Diethylenetriaminepentaacetate",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "tpen": LigandInfo(
        smiles="c1ccc(CN(CCN(Cc2ccccn2)Cc2ccccn2)Cc2ccccn2)nc1",
        mapped_smiles="[cH:1]1[cH:2][cH:3][c:4]([CH2:5][N:6]([CH2:7][CH2:8][N:9]([CH2:10][c:11]2[cH:12][cH:13][cH:14][cH:15][n:16]2)[CH2:17][c:18]2[cH:19][cH:20][cH:21][cH:22][n:23]2)[CH2:24][c:25]2[cH:26][cH:27][cH:28][cH:29][n:30]2)[n:31][cH:32]1",
        denticity=6,
        charge=0,
        aliases=("N,N,N',N'-tetrakis(2-pyridylmethyl)ethylenediamine",),
        description="TPEN",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    # -------------------------------------------------------------------------
    # η-Bonded Ligands
    # -------------------------------------------------------------------------
    "Cp": LigandInfo(
        smiles="c1cc[cH-]c1",
        mapped_smiles="[cH:1]1[cH:2][cH:3][cH-:4][cH:5]1",
        denticity=5,
        charge=-1,
        aliases=("cyclopentadienyl", "C5H5"),
        description="Cyclopentadienyl (η⁵)",
        binding_modes=(
            {
                "name": "eta5",
                "coordination_kind": "haptic",
                "hapticity": 5,
                "donor_mapnums": (1, 2, 3, 4, 5),
                "preferred_bond_type": "dative",
            },
        ),
    ),
    "Cp*": LigandInfo(
        smiles="Cc1c(C)c(C)[c-](C)c1C",
        mapped_smiles="[CH3:1][c:2]1[c:3]([CH3:4])[c:5]([CH3:6])[c-:7]([CH3:8])[c:9]1[CH3:10]",
        denticity=5,
        charge=-1,
        aliases=("pentamethylcyclopentadienyl", "C5Me5", "Cpstar"),
        description="Pentamethylcyclopentadienyl (η⁵)",
        binding_modes=(
            {
                "name": "eta5",
                "coordination_kind": "haptic",
                "hapticity": 5,
                "donor_mapnums": (2, 3, 5, 7, 9),
                "preferred_bond_type": "dative",
            },
        ),
    ),
    # Arenes (η⁶)
    "benzene": LigandInfo(
        smiles="c1ccccc1",
        mapped_smiles="[cH:1]1[cH:2][cH:3][cH:4][cH:5][cH:6]1",
        denticity=6,
        charge=0,
        aliases=("η6-benzene", "C6H6"),
        description="Benzene (η⁶)",
        binding_modes=(
            {
                "name": "eta6",
                "coordination_kind": "haptic",
                "hapticity": 6,
                "donor_mapnums": (1, 2, 3, 4, 5, 6),
                "preferred_bond_type": "dative",
            },
        ),
    ),
    "p-cymene": LigandInfo(
        smiles="Cc1ccc(C(C)C)cc1",
        mapped_smiles="[CH3:1][c:2]1[cH:3][cH:4][c:5]([CH:6]([CH3:7])[CH3:8])[cH:9][cH:10]1",
        denticity=6,
        charge=0,
        aliases=("η6-p-cymene", "4-isopropyltoluene", "p-cym", "p-Cym"),
        description="p-Cymene (η⁶)",
        binding_modes=(
            {
                "name": "eta6",
                "coordination_kind": "haptic",
                "hapticity": 6,
                "donor_mapnums": (2, 3, 4, 5, 9, 10),
                "preferred_bond_type": "dative",
            },
        ),
    ),
    "mesitylene": LigandInfo(
        smiles="Cc1cc(C)cc(C)c1",
        mapped_smiles="[CH3:1][c:2]1[cH:3][c:4]([CH3:5])[cH:6][c:7]([CH3:8])[cH:9]1",
        denticity=6,
        charge=0,
        aliases=("η6-mesitylene", "1,3,5-trimethylbenzene"),
        description="Mesitylene (η⁶)",
        binding_modes=(
            {
                "name": "eta6",
                "coordination_kind": "haptic",
                "hapticity": 6,
                "donor_mapnums": (2, 3, 4, 6, 7, 9),
                "preferred_bond_type": "dative",
            },
        ),
    ),
    "hexamethylbenzene": LigandInfo(
        smiles="Cc1c(C)c(C)c(C)c(C)c1C",
        mapped_smiles="[CH3:1][c:2]1[c:3]([CH3:4])[c:5]([CH3:6])[c:7]([CH3:8])[c:9]([CH3:10])[c:11]1[CH3:12]",
        denticity=6,
        charge=0,
        aliases=("η6-C6Me6", "HMB"),
        description="Hexamethylbenzene (η⁶)",
        binding_modes=(
            {
                "name": "eta6",
                "coordination_kind": "haptic",
                "hapticity": 6,
                "donor_mapnums": (2, 3, 5, 7, 9, 11),
                "preferred_bond_type": "dative",
            },
        ),
    ),
    "toluene": LigandInfo(
        smiles="Cc1ccccc1",
        mapped_smiles="[CH3:1][c:2]1[cH:3][cH:4][cH:5][cH:6][cH:7]1",
        denticity=6,
        charge=0,
        aliases=("η6-toluene",),
        description="Toluene (η⁶)",
        binding_modes=(
            {
                "name": "eta6",
                "coordination_kind": "haptic",
                "hapticity": 6,
                "donor_mapnums": (2, 3, 4, 5, 6, 7),
                "preferred_bond_type": "dative",
            },
        ),
    ),
    "C6H3Me3": LigandInfo(
        smiles="Cc1cc(C)cc(C)c1",
        mapped_smiles="[CH3:1][c:2]1[cH:3][c:4]([CH3:5])[cH:6][c:7]([CH3:8])[cH:9]1",
        denticity=6,
        charge=0,
        aliases=("η6-trimethylbenzene",),
        description="Trimethylbenzene (η⁶)",
        binding_modes=(
            {
                "name": "eta6",
                "coordination_kind": "haptic",
                "hapticity": 6,
                "donor_mapnums": (2, 3, 4, 6, 7, 9),
                "preferred_bond_type": "dative",
            },
        ),
    ),
    # Allyl (η³)
    "η3-allyl": LigandInfo(
        smiles="[CH-]C=C",
        mapped_smiles="[CH-:1][CH:2]=[CH2:3]",
        denticity=3,
        charge=-1,
        aliases=("η3-allyl", "η3-C3H5"),
        description="Allyl (η³)",
        binding_modes=(
            {
                "name": "eta3",
                "coordination_kind": "haptic",
                "hapticity": 3,
                "donor_mapnums": (1, 2, 3),
                "preferred_bond_type": "dative",
            },
        ),
    ),
    "methallyl": LigandInfo(
        smiles="[CH-]C(=C)C",
        mapped_smiles="[CH-:1][C:2](=[CH2:3])[CH3:4]",
        denticity=3,
        charge=-1,
        aliases=("η3-2-methylallyl", "η3-methallyl"),
        description="Methallyl (η³)",
        binding_modes=(
            {
                "name": "eta3",
                "coordination_kind": "haptic",
                "hapticity": 3,
                "donor_mapnums": (1, 2, 3),
                "preferred_bond_type": "dative",
            },
        ),
    ),
    "crotyl": LigandInfo(
        smiles="[CH-]C=CC",
        mapped_smiles="[CH-:1][CH:2]=[CH:3][CH3:4]",
        denticity=3,
        charge=-1,
        aliases=("η3-crotyl", "η3-1-methylallyl"),
        description="Crotyl (η³)",
        binding_modes=(
            {
                "name": "eta3",
                "coordination_kind": "haptic",
                "hapticity": 3,
                "donor_mapnums": (1, 2, 3),
                "preferred_bond_type": "dative",
            },
        ),
    ),
    "cinnamyl": LigandInfo(
        smiles="[CH-]C=Cc1ccccc1",
        mapped_smiles="[CH-:1][CH:2]=[CH:3][c:4]1[cH:5][cH:6][cH:7][cH:8][cH:9]1",
        denticity=3,
        charge=-1,
        aliases=("η3-cinnamyl", "η3-3-phenylallyl"),
        description="Cinnamyl (η³)",
        binding_modes=(
            {
                "name": "eta3",
                "coordination_kind": "haptic",
                "hapticity": 3,
                "donor_mapnums": (1, 2, 3),
                "preferred_bond_type": "dative",
            },
        ),
    ),
    # Tropylium and related
    "C7H7": LigandInfo(
        smiles="[C-]1=CC=CC=CC1",
        mapped_smiles="[C-:1]1=[CH:2][CH:3]=[CH:4][CH:5]=[CH:6][CH2:7]1",
        denticity=7,
        charge=1,
        aliases=("tropylium", "cycloheptatrienyl"),
        description="Tropylium (η⁷)",
        binding_modes=(
            {
                "name": "eta7",
                "coordination_kind": "haptic",
                "hapticity": 7,
                "donor_mapnums": (1, 2, 3, 4, 5, 6, 7),
                "preferred_bond_type": "dative",
            },
        ),
    ),
    "cht": LigandInfo(
        smiles="C1=CC=CCC=C1",
        mapped_smiles="[CH:1]1=[CH:2][CH:3]=[CH:4][CH2:5][CH:6]=[CH:7]1",
        denticity=6,
        charge=0,
        aliases=("cycloheptatriene", "η6-C7H8"),
        description="Cycloheptatriene (η⁶)",
        binding_modes=(
            {
                "name": "eta6",
                "coordination_kind": "haptic",
                "hapticity": 6,
                "donor_mapnums": (1, 2, 3, 4, 6, 7),
                "preferred_bond_type": "dative",
            },
        ),
    ),
    # Indenyl
    "Ind": LigandInfo(
        smiles="[C-]1C=Cc2ccccc21",
        mapped_smiles="[C-:1]1[CH:2]=[CH:3][c:4]2[cH:5][cH:6][cH:7][cH:8][c:9]21",
        denticity=5,
        charge=-1,
        aliases=("indenyl", "η5-indenyl"),
        description="Indenyl (η⁵)",
        binding_modes=(
            {
                "name": "eta5",
                "coordination_kind": "haptic",
                "hapticity": 5,
                "donor_mapnums": (1, 2, 3, 4, 9),
                "preferred_bond_type": "dative",
            },
        ),
    ),
    # Fluorenyl
    "Flu": LigandInfo(
        smiles="[c-]1cccc2c1Cc1ccccc1-2",
        mapped_smiles="[c-:1]1[cH:2][cH:3][cH:4][c:5]2[c:6]1[CH2:7][c:8]1[cH:9][cH:10][cH:11][cH:12][c:13]1-2",
        denticity=5,
        charge=-1,
        aliases=("fluorenyl", "η5-fluorenyl"),
        description="Fluorenyl (η⁵)",
        binding_modes=(
            {
                "name": "eta5",
                "coordination_kind": "haptic",
                "hapticity": 5,
                "donor_mapnums": (5, 6, 7, 8, 13),
                "preferred_bond_type": "dative",
            },
        ),
    ),
    # Cyclooctatetraene
    "cot": LigandInfo(
        smiles="",
        mapped_smiles="",
        denticity=8,
        charge=-2,
        aliases=("cyclooctatetraene", "η8-C8H8", "COT"),
        description="Cyclooctatetraene (η⁸)",
        binding_modes=(
            {
                "name": "eta8",
                "coordination_kind": "haptic",
                "hapticity": 8,
                "donor_mapnums": (1, 2, 3, 4, 5, 6, 7, 8),
                "preferred_bond_type": "dative",
            },
        ),
    ),
    # Butadiene
    "bd": LigandInfo(
        smiles="C=CC=C",
        mapped_smiles="[CH2:1]=[CH:2][CH:3]=[CH2:4]",
        denticity=4,
        charge=0,
        aliases=("butadiene", "η4-butadiene"),
        description="1,3-Butadiene (η⁴)",
        binding_modes=(
            {
                "name": "eta4",
                "coordination_kind": "haptic",
                "hapticity": 4,
                "donor_mapnums": (1, 2, 3, 4),
                "preferred_bond_type": "dative",
            },
        ),
    ),
    "isoprene": LigandInfo(
        smiles="C=CC(=C)C",
        mapped_smiles="[CH2:1]=[CH:2][C:3](=[CH2:4])[CH3:5]",
        denticity=4,
        charge=0,
        aliases=("η4-isoprene", "2-methylbutadiene"),
        description="Isoprene (η⁴)",
        binding_modes=(
            {
                "name": "eta4",
                "coordination_kind": "haptic",
                "hapticity": 4,
                "donor_mapnums": (1, 2, 3, 4),
                "preferred_bond_type": "dative",
            },
        ),
    ),
}


COUNTER_ION_DATABASE: Dict[str, LigandInfo] = {
    "PF6": LigandInfo(
        smiles="F[P-](F)(F)(F)(F)F",
        mapped_smiles="[F:1][P-:2]([F:3])([F:4])([F:5])([F:6])[F:7]",
        charge=-1,
        aliases=("hexafluorophosphate",),
        description="Hexafluorophosphate",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "BF4": LigandInfo(
        smiles="F[B-](F)(F)F",
        mapped_smiles="[F:1][B-:2]([F:3])([F:4])[F:5]",
        charge=-1,
        aliases=("tetrafluoroborate",),
        description="Tetrafluoroborate",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "OTf": LigandInfo(
        smiles="O=S(=O)([O-])C(F)(F)F",
        mapped_smiles="[O:1]=[S:2](=[O:3])([O-:4])[C:5]([F:6])([F:7])[F:8]",
        charge=-1,
        aliases=("triflate", "trifluoromethanesulfonate", "CF3SO3"),
        description="Triflate",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "ClO4": LigandInfo(
        smiles="[O-][Cl+3]([O-])([O-])[O-]",
        mapped_smiles="[O-:1][Cl+3:2]([O-:3])([O-:4])[O-:5]",
        charge=-1,
        aliases=("perchlorate",),
        description="Perchlorate",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "SbF6": LigandInfo(
        smiles="[F][Sb-]([F])([F])([F])([F])[F]",
        mapped_smiles="[F:1][Sb-:2]([F:3])([F:4])([F:5])([F:6])[F:7]",
        charge=-1,
        aliases=("hexafluoroantimonate",),
        description="Hexafluoroantimonate",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "BArF": LigandInfo(
        smiles="",
        mapped_smiles="",
        charge=-1,
        aliases=("BArF24", "tetrakis(3,5-bis(trifluoromethyl)phenyl)borate"),
        description="Tetrakis(3,5-bis(trifluoromethyl)phenyl)borate",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "BAr4": LigandInfo(
        smiles="c1ccc([B-](c2ccccc2)(c2ccccc2)c2ccccc2)cc1",
        mapped_smiles="[cH:1]1[cH:2][cH:3][c:4]([B-:5]([c:6]2[cH:7][cH:8][cH:9][cH:10][cH:11]2)([c:12]2[cH:13][cH:14][cH:15][cH:16][cH:17]2)[c:18]2[cH:19][cH:20][cH:21][cH:22][cH:23]2)[cH:24][cH:25]1",
        charge=-1,
        aliases=("tetraphenylborate", "BPh4"),
        description="Tetraphenylborate",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "NO3": LigandInfo(
        smiles="O=[N+]([O-])[O-]",
        mapped_smiles="[O:1]=[N+:2]([O-:3])[O-:4]",
        charge=-1,
        aliases=("nitrate",),
        description="Nitrate",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "Cl": LigandInfo(
        smiles="[Cl-]",
        mapped_smiles="[Cl-:1]",
        charge=-1,
        aliases=("chloride",),
        description="Chloride",
        binding_modes=(
            {
                "name": "kappa-Cl",
                "coordination_kind": "kappa",
                "hapticity": None,
                "donor_mapnums": (1,),
                "preferred_bond_type": "single",
            },
        ),
    ),
    "Br": LigandInfo(
        smiles="[Br-]",
        mapped_smiles="[Br-:1]",
        charge=-1,
        aliases=("bromide",),
        description="Bromide",
        binding_modes=(
            {
                "name": "kappa-Br",
                "coordination_kind": "kappa",
                "hapticity": None,
                "donor_mapnums": (1,),
                "preferred_bond_type": "single",
            },
        ),
    ),
    "I": LigandInfo(
        smiles="[I-]",
        mapped_smiles="[I-:1]",
        charge=-1,
        aliases=("iodide",),
        description="Iodide",
        binding_modes=(
            {
                "name": "kappa-I",
                "coordination_kind": "kappa",
                "hapticity": None,
                "donor_mapnums": (1,),
                "preferred_bond_type": "single",
            },
        ),
    ),
    "BPh4": LigandInfo(
        smiles="c1ccc([B-](c2ccccc2)(c2ccccc2)c2ccccc2)cc1",
        mapped_smiles="[cH:1]1[cH:2][cH:3][c:4]([B-:5]([c:6]2[cH:7][cH:8][cH:9][cH:10][cH:11]2)([c:12]2[cH:13][cH:14][cH:15][cH:16][cH:17]2)[c:18]2[cH:19][cH:20][cH:21][cH:22][cH:23]2)[cH:24][cH:25]1",
        charge=-1,
        aliases=("tetraphenylborate",),
        description="Tetraphenylborate",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "Al(OC(CF3)3)4": LigandInfo(
        smiles="",
        mapped_smiles="",
        charge=-1,
        aliases=("perfluoro-tert-butoxide aluminate",),
        description="Perfluoro-tert-butoxide aluminate",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "BAr4F": LigandInfo(
        smiles="",
        mapped_smiles="",
        charge=-1,
        aliases=("tetrakis(3,5-bis(trifluoromethyl)phenyl)borate", "BArF"),
        description="BArF",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "NTf2": LigandInfo(
        smiles="O=S(=O)([N-]S(=O)(=O)C(F)(F)F)C(F)(F)F",
        mapped_smiles="[O:1]=[S:2](=[O:3])([N-:4][S:5](=[O:6])(=[O:7])[C:8]([F:9])([F:10])[F:11])[C:12]([F:13])([F:14])[F:15]",
        charge=-1,
        aliases=("bis(trifluoromethylsulfonyl)imide", "TFSI", "bistriflimide"),
        description="Bis(trifluoromethylsulfonyl)imide",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "AsF6": LigandInfo(
        smiles="F[As-](F)(F)(F)(F)F",
        mapped_smiles="[F:1][As-:2]([F:3])([F:4])([F:5])([F:6])[F:7]",
        charge=-1,
        aliases=("hexafluoroarsenate",),
        description="Hexafluoroarsenate",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "B(C6F5)4": LigandInfo(
        smiles="Fc1c(F)c(F)c([B-](c2c(F)c(F)c(F)c(F)c2F)(c2c(F)c(F)c(F)c(F)c2F)c2c(F)c(F)c(F)c(F)c2F)c(F)c1F",
        mapped_smiles="[F:1][c:2]1[c:3]([F:4])[c:5]([F:6])[c:7]([B-:8]([c:9]2[c:10]([F:11])[c:12]([F:13])[c:14]([F:15])[c:16]([F:17])[c:18]2[F:19])([c:20]2[c:21]([F:22])[c:23]([F:24])[c:25]([F:26])[c:27]([F:28])[c:29]2[F:30])[c:31]2[c:32]([F:33])[c:34]([F:35])[c:36]([F:37])[c:38]([F:39])[c:40]2[F:41])[c:42]([F:43])[c:44]1[F:45]",
        charge=-1,
        aliases=("tetrakis(pentafluorophenyl)borate",),
        description="Tetrakis(pentafluorophenyl)borate",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "HSO4": LigandInfo(
        smiles="O=S(=O)([O-])O",
        mapped_smiles="[O:1]=[S:2](=[O:3])([O-:4])[OH:5]",
        charge=-1,
        aliases=("hydrogensulfate", "bisulfate"),
        description="Hydrogen sulfate",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "CF3CO2": LigandInfo(
        smiles="O=C([O-])C(F)(F)F",
        mapped_smiles="[O:1]=[C:2]([O-:3])[C:4]([F:5])([F:6])[F:7]",
        charge=-1,
        aliases=("trifluoroacetate", "TFA"),
        description="Trifluoroacetate",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "MeSO3": LigandInfo(
        smiles="CS(=O)(=O)[O-]",
        mapped_smiles="[CH3:1][S:2](=[O:3])(=[O:4])[O-:5]",
        charge=-1,
        aliases=("mesylate", "methanesulfonate"),
        description="Mesylate",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "TsO": LigandInfo(
        smiles="Cc1ccc(S(=O)(=O)[O-])cc1",
        mapped_smiles="[CH3:1][c:2]1[cH:3][cH:4][c:5]([S:6](=[O:7])(=[O:8])[O-:9])[cH:10][cH:11]1",
        charge=-1,
        aliases=("tosylate", "4-toluenesulfonate", "OTs"),
        description="Tosylate",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "ReO4": LigandInfo(
        smiles="[O]=[Re]([OH])([OH])([OH])[OH]",
        mapped_smiles="[O:1]=[Re:2]([OH:3])([OH:4])([OH:5])[OH:6]",
        charge=-1,
        aliases=("perrhenate",),
        description="Perrhenate",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "IO4": LigandInfo(
        smiles="OI(O)(O)(O)O",
        mapped_smiles="[OH:1][I:2]([OH:3])([OH:4])([OH:5])[OH:6]",
        charge=-1,
        aliases=("periodate",),
        description="Periodate",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "BH4": LigandInfo(
        smiles="[BH4-]",
        mapped_smiles="[BH4-:1]",
        charge=-1,
        aliases=("borohydride", "tetrahydroborate"),
        description="Borohydride (as counter ion)",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "AlH4": LigandInfo(
        smiles="[AlH4-]",
        mapped_smiles="[AlH4-:1]",
        charge=-1,
        aliases=("aluminate", "tetrahydroaluminate"),
        description="Tetrahydroaluminate",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "AlCl4": LigandInfo(
        smiles="[Cl][Al-]([Cl])([Cl])[Cl]",
        mapped_smiles="[Cl:1][Al-:2]([Cl:3])([Cl:4])[Cl:5]",
        charge=-1,
        aliases=("tetrachloroaluminate",),
        description="Tetrachloroaluminate",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "FeCl4": LigandInfo(
        smiles="[Cl][Fe-]([Cl])([Cl])[Cl]",
        mapped_smiles="[Cl:1][Fe-:2]([Cl:3])([Cl:4])[Cl:5]",
        charge=-1,
        aliases=("tetrachloroferrate",),
        description="Tetrachloroferrate",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "CuCl2": LigandInfo(
        smiles="[Cl][Cu-][Cl]",
        mapped_smiles="[Cl:1][Cu-:2][Cl:3]",
        charge=-1,
        aliases=("dichlorocuprate",),
        description="Dichlorocuprate",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "ZnCl3": LigandInfo(
        smiles="[Cl][Zn-]([Cl])[Cl]",
        mapped_smiles="[Cl:1][Zn-:2]([Cl:3])[Cl:4]",
        charge=-1,
        aliases=("trichlorozincate",),
        description="Trichlorozincate",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "GaCl4": LigandInfo(
        smiles="[Cl][Ga-]([Cl])([Cl])[Cl]",
        mapped_smiles="[Cl:1][Ga-:2]([Cl:3])([Cl:4])[Cl:5]",
        charge=-1,
        aliases=("tetrachlorogallate",),
        description="Tetrachlorogallate",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    # Cationic counterions (for anionic complexes)
    "Na": LigandInfo(
        smiles="[Na+]",
        mapped_smiles="[Na+:1]",
        charge=1,
        aliases=("sodium",),
        description="Sodium",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "K": LigandInfo(
        smiles="[K+]",
        mapped_smiles="[K+:1]",
        charge=1,
        aliases=("potassium",),
        description="Potassium",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "Li": LigandInfo(
        smiles="[Li+]",
        mapped_smiles="[Li+:1]",
        charge=1,
        aliases=("lithium",),
        description="Lithium",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "Cs": LigandInfo(
        smiles="[Cs+]",
        mapped_smiles="[Cs+:1]",
        charge=1,
        aliases=("cesium", "caesium"),
        description="Cesium",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "NBu4": LigandInfo(
        smiles="CCCC[N+](CCCC)(CCCC)CCCC",
        mapped_smiles="[CH3:1][CH2:2][CH2:3][CH2:4][N+:5]([CH2:6][CH2:7][CH2:8][CH3:9])([CH2:10][CH2:11][CH2:12][CH3:13])[CH2:14][CH2:15][CH2:16][CH3:17]",
        charge=1,
        aliases=("tetrabutylammonium", "TBA", "nBu4N"),
        description="Tetrabutylammonium",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "NEt4": LigandInfo(
        smiles="CC[N+](CC)(CC)CC",
        mapped_smiles="[CH3:1][CH2:2][N+:3]([CH2:4][CH3:5])([CH2:6][CH3:7])[CH2:8][CH3:9]",
        charge=1,
        aliases=("tetraethylammonium", "TEA", "Et4N"),
        description="Tetraethylammonium",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "NMe4": LigandInfo(
        smiles="C[N+](C)(C)C",
        mapped_smiles="[CH3:1][N+:2]([CH3:3])([CH3:4])[CH3:5]",
        charge=1,
        aliases=("tetramethylammonium", "TMA", "Me4N"),
        description="Tetramethylammonium",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "PPh4": LigandInfo(
        smiles="c1ccc([P+](c2ccccc2)(c2ccccc2)c2ccccc2)cc1",
        mapped_smiles="[cH:1]1[cH:2][cH:3][c:4]([P+:5]([c:6]2[cH:7][cH:8][cH:9][cH:10][cH:11]2)([c:12]2[cH:13][cH:14][cH:15][cH:16][cH:17]2)[c:18]2[cH:19][cH:20][cH:21][cH:22][cH:23]2)[cH:24][cH:25]1",
        charge=1,
        aliases=("tetraphenylphosphonium", "Ph4P"),
        description="Tetraphenylphosphonium",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "PPN": LigandInfo(
        smiles="",
        mapped_smiles="",
        charge=1,
        aliases=("bis(triphenylphosphine)iminium", "Ph3P=N=PPh3"),
        description="Bis(triphenylphosphine)iminium",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "Cp2Fe": LigandInfo(
        smiles="",
        mapped_smiles="",
        charge=1,
        aliases=("ferrocenium", "Fc+"),
        description="Ferrocenium",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "Cp2Co": LigandInfo(
        smiles="C1=CCC=C1.[Co+2].c1cc[cH-]c1",
        mapped_smiles="[CH:1]1=[CH:2][CH2:3][CH:4]=[CH:5]1.[Co+2:6].[cH:7]1[cH:8][cH:9][cH-:10][cH:11]1",
        charge=1,
        aliases=("cobaltocenium",),
        description="Cobaltocenium",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "H3O": LigandInfo(
        smiles="",
        mapped_smiles="",
        charge=1,
        aliases=("hydronium", "oxonium"),
        description="Hydronium",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "NH4": LigandInfo(
        smiles="[NH4+]",
        mapped_smiles="[NH4+:1]",
        charge=1,
        aliases=("ammonium",),
        description="Ammonium",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "pyH": LigandInfo(
        smiles="c1cc[nH+]cc1",
        mapped_smiles="[cH:1]1[cH:2][cH:3][nH+:4][cH:5][cH:6]1",
        charge=1,
        aliases=("pyridinium",),
        description="Pyridinium",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
    "DMAH": LigandInfo(
        smiles="C[NH+](C)c1ccccc1",
        mapped_smiles="[CH3:1][NH+:2]([CH3:3])[c:4]1[cH:5][cH:6][cH:7][cH:8][cH:9]1",
        charge=1,
        aliases=("dimethylanilinium",),
        description="Dimethylanilinium",
        binding_modes=(
            {
                "name": "",
                "coordination_kind": "",
                "hapticity": None,
                "donor_mapnums": (),
                "preferred_bond_type": "",
            },
        ),
    ),
}


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
