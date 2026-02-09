from typing import Dict, List

import pubchempy as pcp

from cholla_chem.utils.logging_config import logger
from cholla_chem.utils.string_utils import filter_latin1_compatible


def name_to_smiles_pubchem(compound_name_list: List[str]) -> Dict[str, str]:
    """
    Convert chemical names to SMILES using pubchem.

    Args:
        compound_name_list (List[str]): List of compound names to convert to SMILES.

    Returns:
        Dict[str, str]: Dictionary of compound names to SMILES.
    """
    # Filter out non-latin1 compatible strings
    filtered_compound_name_list = filter_latin1_compatible(compound_name_list)
    if not filtered_compound_name_list:
        return {}

    pubchem_compounds = pcp.get_compounds(filtered_compound_name_list, "name")
    pubchem_name_dict = {
        k: v.smiles if v is not None else ""
        for k, v in zip(filtered_compound_name_list, pubchem_compounds)
    }

    if len(filtered_compound_name_list) != len(pubchem_name_dict):
        logger.warning(
            f"Mismatching lengths: "
            f"filtered_compound_name_list ({len(filtered_compound_name_list)}), "
            f"pubchem_name_dict ({len(pubchem_name_dict)})"
        )
        print(filtered_compound_name_list)
        print([ele.smiles if ele else "" for ele in pubchem_compounds])
        return {}

    return pubchem_name_dict
