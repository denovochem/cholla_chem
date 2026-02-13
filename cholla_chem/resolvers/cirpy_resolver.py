from typing import Dict, List

import cirpy

from cholla_chem.utils.logging_config import logger


def retrieve_cirpy_results(compound_name: str) -> str:
    """
    Retrieves the SMILES string for a given compound identifier using CIRpy.

    Args:
        compound_name (str): The compound name.

    Returns:
        str: The SMILES string of the compound, or an empty string if an exception occurs.
    """
    try:
        smiles = cirpy.resolve(compound_name, "smiles")

    except Exception as e:
        logger.warning(f"Exception with CIRpy query: {str(e)}")
        return ""

    return smiles


def name_to_smiles_cirpy(compound_name_list: List[str]) -> Dict[str, str]:
    """
    Converts a list of chemical names to their corresponding SMILES strings using CIRpy.

    Args:
        compound_name_list (List[str]): A list of chemical names to be converted.

    Returns:
        Dict[str, str]: A dictionary mapping each chemical name to its SMILES string.
    """
    cirpy_name_dict = {}
    for compound_name in compound_name_list:
        result = retrieve_cirpy_results(compound_name)
        if result:
            cirpy_name_dict[compound_name] = result
    return cirpy_name_dict
