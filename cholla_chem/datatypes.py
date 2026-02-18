from typing import Dict, List, Tuple, TypedDict


class CompoundResolutionEntry(TypedDict):
    SMILES: str
    SMILES_source: List[str]
    SMILES_dict: Dict[str, List[str]]  # canonical_smiles -> [resolver, ...]
    additional_info: Dict[str, str]  # resolver -> info message


NameCorrectionCandidate = Tuple[str, float]


class NameCorrectionInfo(TypedDict):
    top_5: List[NameCorrectionCandidate]
    selected_name: str
    name_manipulation_method: str
    SMILES: str


class CompoundResolutionEntryWithNameCorrection(TypedDict):
    SMILES: str
    SMILES_source: List[str]
    SMILES_dict: Dict[str, List[str]]  # canonical_smiles -> [resolver, ...]
    additional_info: Dict[str, str]  # resolver -> info message
    name_correction_info: NameCorrectionInfo
