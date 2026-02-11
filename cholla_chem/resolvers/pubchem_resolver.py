from typing import Dict, List

import pubchempy as pcp
import requests

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

    try:
        pubchem_compounds = pcp.get_compounds(filtered_compound_name_list, "name")
    except requests.exceptions.HTTPError as http_err:
        logger.warning(f"HTTP error in PubChem query: {http_err}")
        return {}
    except requests.exceptions.RequestException as err:
        logger.warning(f"Request error in PubChem query: {err}")
        return {}
    except Exception as e:
        logger.warning(f"Unexpected error in PubChem query: {e}")
        return {}

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
        return {}

    return pubchem_name_dict
