from __future__ import annotations

import json
import os
import ssl
from typing import Any, Dict, List
from urllib.error import HTTPError, URLError
from urllib.parse import quote
from urllib.request import urlopen

from cholla_chem.utils.logging_config import logger
from cholla_chem.utils.string_utils import filter_latin1_compatible

_CA_FILE = os.getenv("PUBCHEMPY_CA_BUNDLE") or os.getenv("REQUESTS_CA_BUNDLE")
if not _CA_FILE:
    try:
        import certifi

        _CA_FILE = certifi.where()
    except ImportError:
        _CA_FILE = None

API_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"


def request_smiles_by_name(compound_name: str, timeout: int = 15) -> bytes:
    """
    Request the canonical SMILES for a single compound name from PubChem.

    Args:
        compound_name (str): Compound name to resolve.
        timeout (int): Timeout in seconds for the HTTP request.

    Returns:
        bytes: Raw response body from the PubChem PUG REST API.

    Raises:
        ValueError: If `compound_name` is empty.
        Exception: If the PubChem request fails.
    """
    if not compound_name:
        raise ValueError("compound_name cannot be empty")

    encoded_name = quote(compound_name, safe="")
    api_url = f"{API_BASE}/compound/name/{encoded_name}/property/SMILES/JSON"
    context = ssl.create_default_context(cafile=_CA_FILE)

    try:
        with urlopen(api_url, timeout=timeout, context=context) as response:
            return response.read()
    except HTTPError as e:
        raise Exception(f"HTTP Error {e.code}: {e.msg}") from e
    except URLError as e:
        raise Exception(f"URL Error: {e.reason}") from e


def get_smiles_by_name(compound_name: str) -> str:
    """
    Resolve a single compound name to a canonical SMILES string.

    Args:
        compound_name (str): Compound name to resolve.

    Returns:
        str: The resolved canonical SMILES string, or an empty string if resolution
            fails or no SMILES is available.
    """
    try:
        response = request_smiles_by_name(compound_name)
        payload: Dict[str, Any] = json.loads(response.decode())
    except Exception as e:
        logger.info(e)
        return ""

    property_table = payload.get("PropertyTable", {})
    properties = property_table.get("Properties", [])
    if not properties:
        return ""

    first_result = properties[0]
    smiles = first_result.get("SMILES", "")
    if not isinstance(smiles, str):
        return ""

    return smiles


def name_to_smiles_pubchem(compound_name_list: List[str]) -> Dict[str, str]:
    """
    Convert chemical names to SMILES using PubChem PUG REST single-name lookups.

    Args:
        compound_name_list (List[str]): List of compound names to convert to SMILES.

    Returns:
        Dict[str, str]: Dictionary of compound names to SMILES.
    """
    filtered_compound_name_list = filter_latin1_compatible(compound_name_list)
    if not filtered_compound_name_list:
        return {}

    pubchem_name_dict: Dict[str, str] = {}
    for compound_name in filtered_compound_name_list:
        pubchem_name_dict[compound_name] = get_smiles_by_name(compound_name)

    if len(filtered_compound_name_list) != len(pubchem_name_dict):
        logger.warning(
            f"Mismatching lengths: "
            f"filtered_compound_name_list ({len(filtered_compound_name_list)}), "
            f"pubchem_name_dict ({len(pubchem_name_dict)})"
        )
        return {}

    return pubchem_name_dict
