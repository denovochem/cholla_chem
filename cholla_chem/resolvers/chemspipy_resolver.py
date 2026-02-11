from typing import Dict, List

import requests
from chemspipy import ChemSpider

from cholla_chem.utils.logging_config import logger
from cholla_chem.utils.string_utils import filter_latin1_compatible


def name_to_smiles_chemspipy(
    compound_name_list: List[str],
    chemspider_api_key: str,
) -> Dict[str, str]:
    """
    Convert chemical names to SMILES using ChemSpiPy.

    Args:
        compound_name_list (List[str]): List of compound names to convert to SMILES.
        chemspider_api_key (str): ChemSpider API key (https://developer.rsc.org/getting-started)

    Returns:
        Dict[str, str]: Dictionary of compound names to SMILES.
    """
    try:
        cs = ChemSpider(chemspider_api_key)
    except Exception as e:
        logger.warning(f"Error initializing ChemSpiPy: {e}")
        return {}

    chemspipy_name_dict = {}
    for compound_name in filter_latin1_compatible(compound_name_list):
        try:
            results = cs.search(compound_name)
            results.wait()
            results.ready()
            if len(results) == 0:
                continue
            c = results[0]
            if not c.smiles:
                continue
            chemspipy_name_dict[compound_name] = c.smiles
        except requests.exceptions.HTTPError as http_err:
            logger.warning(f"HTTP error in ChemSpiPy query: {http_err}")
            continue
        except requests.exceptions.RequestException as err:
            logger.warning(f"Request error in ChemSpiPy query: {err}")
            continue

    return chemspipy_name_dict
