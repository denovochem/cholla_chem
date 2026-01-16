from typing import Dict, List

from placeholder_name.name_manipulation.name_correction.name_corrector import (
    ChemNameCorrector,
)
from placeholder_name.name_manipulation.peptide_shorthand_handler import (
    peptide_shorthand_to_iupac,
    looks_like_peptide_shorthand,
)


corrector = ChemNameCorrector()


def correct_names(names_to_correct: list[str]) -> dict[str, Dict[str, List[str]]]:
    """
    Corrects a list of chemical names using OCR and Typo correction strategies.

    Args:
        names_to_correct (list[str]): List of chemical names to correct

    Returns:
        dict[str, Dict[str,List[str]]]: A dictionary mapping each original name to its expanded and/or corrected names. The value is a dictionary with keys "expanded" and/or "corrected", and values being dictionaries with the expanded/corrected name as key and a list of operations used to derive it as value.
    """
    original_name_corrected_name_dict: dict[str, Dict[str, List[str]]] = {}
    if not names_to_correct:
        return original_name_corrected_name_dict
    for name_to_correct in names_to_correct:
        if looks_like_peptide_shorthand(name_to_correct):
            if name_to_correct not in original_name_corrected_name_dict:
                original_name_corrected_name_dict[name_to_correct] = {
                    peptide_shorthand_to_iupac(name_to_correct): [
                        "peptide_shorthand_expansion"
                    ]
                }

        if corrector.correct(name_to_correct) is not None:
            if name_to_correct not in original_name_corrected_name_dict:
                original_name_corrected_name_dict[name_to_correct] = {
                    corrector.correct(name_to_correct): ["name_correction"]
                }

    return original_name_corrected_name_dict
