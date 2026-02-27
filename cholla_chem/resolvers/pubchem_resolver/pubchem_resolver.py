from __future__ import annotations

import json
import os
import ssl
import time
import xml.etree.ElementTree as ET
from http.client import HTTPResponse
from typing import Any, Dict, List, Optional
from urllib.error import HTTPError
from urllib.parse import quote, urlencode
from urllib.request import Request, urlopen

from cholla_chem.utils.logging_config import logger
from cholla_chem.utils.string_utils import filter_latin1_compatible

# Get SSL certs from env var or certifi package if available.
_CA_FILE = os.getenv("PUBCHEMPY_CA_BUNDLE") or os.getenv("REQUESTS_CA_BUNDLE")
if not _CA_FILE:
    try:
        import certifi

        _CA_FILE = certifi.where()
    except ImportError:
        _CA_FILE = None

# Dummy values for instances where pubchem fails to return a result
# Selected to be a valid but exceedingly rare CID and SMILES
DUMMY_CID = "58965162"
DUMMY_SMILES = "[22CH4]"

# Base URL for the PubChem PUG REST API.
API_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"

# URL for batch identifier exchange service
PUBCHEM_PUG_URL = "https://pubchem.ncbi.nlm.nih.gov/pug/pug.cgi"

# Type alias for URL query parameters.
QueryParam = str | int | float | bool | list[str] | None


def request(
    identifier: str | int | List[str] | List[int] | List[str | int],
    namespace: str = "cid",
    domain: str = "compound",
    operation: str | None = None,
    output: str = "JSON",
    searchtype: str | None = None,
    **kwargs: QueryParam,
) -> HTTPResponse:
    """Construct API request from parameters and return the response.

    Full specification at https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest
    """
    if not identifier:
        raise ValueError("identifier/cid cannot be None")
    # If identifier is an integer, convert to string
    if isinstance(identifier, int):
        identifier = str(identifier)
    # If identifier is a list and in specified namespace, batch retrieve cids, then join with commas into string
    if not isinstance(identifier, str) and namespace in [
        "name",
        "smiles",
        "inchi",
        "inchikey",
    ]:
        id_list = [str(x) for x in identifier]
        cid_list = batch_retrieve_cids(id_list, namespace)
        identifier = ",".join(cid_list)
        namespace = "cid"
    # If identifier is a list, join with commas into string
    elif not isinstance(identifier, str):
        identifier = ",".join(str(x) for x in identifier)
    # Filter None values from kwargs
    kwargs = {k: v for k, v in kwargs.items() if v is not None}
    # Build API URL
    urlid, postdata = None, None
    if namespace == "sourceid":
        identifier = identifier.replace("/", ".")
    if (
        namespace in ["listkey", "formula", "sourceid"]
        or searchtype == "xref"
        or (searchtype and namespace == "cid")
        or domain == "sources"
    ):
        urlid = quote(identifier.encode("utf8"))
    else:
        postdata = urlencode([(namespace, identifier)]).encode("utf8")
    comps = filter(
        None, [API_BASE, domain, searchtype, namespace, urlid, operation, output]
    )
    apiurl = "/".join(comps)
    if kwargs:
        apiurl += f"?{urlencode(kwargs)}"
    # Make request
    try:
        context = ssl.create_default_context(cafile=_CA_FILE)
        response = urlopen(apiurl, postdata, context=context)
        return response
    except HTTPError as e:
        raise Exception(f"HTTP Error {e.code}: {e.msg}") from e


def get(
    identifier: str | int | List[str] | List[int] | List[str | int],
    namespace: str = "cid",
    domain: str = "compound",
    operation: str | None = None,
    output: str = "JSON",
    searchtype: str | None = None,
    **kwargs: QueryParam,
) -> bytes:
    """Request wrapper that automatically handles async requests."""
    if (searchtype and searchtype != "xref") or namespace in ["formula"]:
        response = request(
            identifier, namespace, domain, None, "JSON", searchtype, **kwargs
        ).read()
        status = json.loads(response.decode())
        if "Waiting" in status and "ListKey" in status["Waiting"]:
            identifier = status["Waiting"]["ListKey"]
            namespace = "listkey"
            while "Waiting" in status and "ListKey" in status["Waiting"]:
                time.sleep(2)
                response = request(
                    identifier,
                    namespace=namespace,
                    domain=domain,
                    operation=operation,
                    output="JSON",
                    searchtype=None,
                    **kwargs,
                ).read()

                status = json.loads(response.decode())
            if not output == "JSON":
                response = request(
                    identifier,
                    namespace,
                    domain,
                    operation,
                    output,
                    searchtype,
                    **kwargs,
                ).read()
    else:
        response = request(
            identifier, namespace, domain, operation, output, searchtype, **kwargs
        ).read()
    return response


def batch_retrieve_cids(
    identifiers: List[str],
    namespace: str = "name",
    chunk_size: int = 1000,
    check_interval: int = 3,
    timeout: int = 120,
) -> List[str]:
    """Retrieves cids from PubChem in batches.

    This is the main function for batch retrieval. It chunks the input identifiers,
    submits batch requests, polls for completion, and returns aggregated results.
    """

    # Validate chunk size
    if chunk_size > 1000:
        chunk_size = 1000
        logger.warning("Chunk size limited to 1000")

    # Process in chunks
    cids = []
    for i in range(0, len(identifiers), chunk_size):
        chunk = identifiers[i : i + chunk_size]

        # Create and send the batch request
        xml_request = create_batch_cid_request_xml(chunk, namespace)
        initial_response = send_xml_to_pug(xml_request)

        # Extract request ID
        root = ET.fromstring(initial_response)
        req_id_elem = root.find(".//PCT-Waiting_reqid")

        if req_id_elem is None:
            # Check for immediate error in response
            error_elem = root.find(".//PCT-Status[@value='error']")
            if error_elem is not None:
                error_message_elem = root.find(".//PCT-Status-Message_message")
                error_message = (
                    error_message_elem.text
                    if error_message_elem is not None
                    and error_message_elem.text is not None
                    else "Unknown error"
                )
                logger.warning(f"Error retrieving CIDs from PubChem: {error_message}")

            cids.extend([DUMMY_CID] * len(chunk))
            continue

        req_id = req_id_elem.text
        if not req_id:
            cids.extend([DUMMY_CID] * len(chunk))
            continue

        # Poll for completion
        download_url = poll_request_status(req_id, check_interval, timeout)

        if download_url is None:
            cids.extend([DUMMY_CID] * len(chunk))
            continue

        # Download and parse CID file
        file_content = download_file_from_pug(download_url)
        identifier_to_cid = parse_cid_file(file_content)

        if not identifier_to_cid:
            cids.extend([DUMMY_CID] * len(chunk))
            continue

        # Fetch CIDs
        chunk_cids = [
            identifier_to_cid.get(identifier, DUMMY_CID) for identifier in chunk
        ]
        cids.extend(chunk_cids)

    return cids


def send_xml_to_pug(xml_data: str) -> str:
    """Sends an XML request to the PubChem PUG REST service."""
    headers = {"Content-Type": "application/xml"}
    data = xml_data.encode("utf-8")
    req = Request(PUBCHEM_PUG_URL, data=data, headers=headers, method="POST")
    context = ssl.create_default_context(cafile=_CA_FILE)

    try:
        with urlopen(req, timeout=30, context=context) as resp:
            body = resp.read()
            return body.decode("utf-8", errors="replace")
    except HTTPError as e:
        raise Exception(f"HTTP Error {e.code}: {e.msg}") from e


def create_batch_cid_request_xml(
    identifiers: List[str],
    identifier_type: str,
) -> str:
    """Creates the initial XML request for batch property retrieval using ID exchange."""
    root = ET.Element("PCT-Data")
    input_data = ET.SubElement(root, "PCT-Data_input")
    input_data = ET.SubElement(input_data, "PCT-InputData")

    # Query section
    query = ET.SubElement(input_data, "PCT-InputData_query")
    query = ET.SubElement(query, "PCT-Query")
    query_type = ET.SubElement(query, "PCT-Query_type")
    query_type = ET.SubElement(query_type, "PCT-QueryType")

    # ID exchange section
    id_exchange = ET.SubElement(query_type, "PCT-QueryType_id-exchange")
    id_exchange = ET.SubElement(id_exchange, "PCT-QueryIDExchange")

    # Input section with identifiers
    input_elem = ET.SubElement(id_exchange, "PCT-QueryIDExchange_input")
    query_uids = ET.SubElement(input_elem, "PCT-QueryUids")

    # Add identifiers based on type
    if identifier_type.lower() == "name":
        name_elem = ET.SubElement(query_uids, "PCT-QueryUids_synonyms")
        for identifier in identifiers:
            ET.SubElement(name_elem, "PCT-QueryUids_synonyms_E").text = str(identifier)
    elif identifier_type.lower() == "smiles":
        smiles_elem = ET.SubElement(query_uids, "PCT-QueryUids_smiles")
        for identifier in identifiers:
            ET.SubElement(smiles_elem, "PCT-QueryUids_smiles_E").text = str(identifier)
    elif identifier_type.lower() == "inchi":
        inchi_elem = ET.SubElement(query_uids, "PCT-QueryUids_inchis")
        for identifier in identifiers:
            ET.SubElement(inchi_elem, "PCT-QueryUids_inchis_E").text = str(identifier)
    elif identifier_type.lower() == "inchikey":
        inchikey_elem = ET.SubElement(query_uids, "PCT-QueryUids_inchi-keys")
        for identifier in identifiers:
            ET.SubElement(inchikey_elem, "PCT-QueryUids_inchi-keys_E").text = str(
                identifier
            )
    else:
        # Default to synonyms for other types
        synonyms = ET.SubElement(query_uids, "PCT-QueryUids_synonyms")
        for identifier in identifiers:
            ET.SubElement(synonyms, "PCT-QueryUids_synonyms_E").text = str(identifier)

    # Operation and output settings
    ET.SubElement(id_exchange, "PCT-QueryIDExchange_operation-type").set(
        "value", "same"
    )
    ET.SubElement(id_exchange, "PCT-QueryIDExchange_output-type").set("value", "cid")
    ET.SubElement(id_exchange, "PCT-QueryIDExchange_output-method").set(
        "value", "file-pair"
    )
    ET.SubElement(id_exchange, "PCT-QueryIDExchange_compression").set("value", "none")

    return ET.tostring(root, encoding="unicode")


def create_status_request_xml(req_id: str) -> str:
    """Creates the XML request to check the status of a previously submitted request."""
    root = ET.Element("PCT-Data")
    input_data = ET.SubElement(root, "PCT-Data_input")
    input_data = ET.SubElement(input_data, "PCT-InputData")
    request = ET.SubElement(input_data, "PCT-InputData_request")
    request = ET.SubElement(request, "PCT-Request")
    ET.SubElement(request, "PCT-Request_reqid").text = req_id
    ET.SubElement(request, "PCT-Request_type").set("value", "status")

    return ET.tostring(root, encoding="unicode")


def download_file_from_pug(url: str) -> str:
    """Downloads a file from the given URL."""
    # Convert FTP to HTTPS if needed
    https_url = url.replace("ftp://", "https://")
    req = Request(https_url, method="GET")
    context = ssl.create_default_context(cafile=_CA_FILE)

    try:
        with urlopen(req, timeout=60, context=context) as resp:
            body = resp.read()
            return body.decode("utf-8", errors="replace")
    except HTTPError as e:
        raise Exception(f"HTTP Error {e.code}: {e.msg}") from e


def poll_request_status(
    req_id: str, check_interval: int = 3, timeout: int = 120
) -> Optional[str]:
    """Polls PubChem for request completion and returns the download URL."""
    start_time = time.time()

    while time.time() - start_time < timeout:
        status_xml = create_status_request_xml(req_id)

        try:
            status_response = send_xml_to_pug(status_xml)
        except Exception:
            return None

        root = ET.fromstring(status_response)
        status_elem = root.find(".//PCT-Status")

        if status_elem is None:
            return None

        status = status_elem.attrib.get("value", "unknown")

        if status == "success":
            download_url_elem = root.find(".//PCT-Download-URL_url")
            if download_url_elem is not None:
                return download_url_elem.text
            else:
                return None

        elif status == "error":
            error_elem = root.find(".//PCT-Status-Message_message")
            error_message = (
                error_elem.text if error_elem is not None else "Unknown error"
            )
            logger.warning(f"Error retrieving CIDs from PubChem: {error_message}")
            return None

        elif status in ["queued", "running"]:
            time.sleep(check_interval)
        else:
            time.sleep(check_interval)

    return None


def parse_cid_file(file_content: str) -> Dict[str, str]:
    """Parses the CID file from PubChem into a dictionary mapping input to CID."""
    identifier_to_cid = {}
    lines = file_content.strip().split("\n")

    for i, line in enumerate(lines):
        parts = line.split("\t")
        if len(parts) >= 2:
            input_identifier = parts[0].strip()
            cid = parts[1].strip()
            # Only add if CID is numeric (valid)
            if cid and cid.isdigit():
                identifier_to_cid[input_identifier] = cid
            else:
                identifier_to_cid[input_identifier] = DUMMY_CID
        else:
            input_identifier = parts[0].strip()
            identifier_to_cid[input_identifier] = DUMMY_CID

    return identifier_to_cid


def get_json(
    identifier: str | int | List[str] | List[int] | List[str | int],
    namespace: str = "cid",
    domain: str = "compound",
    operation: str | None = None,
    searchtype: str | None = None,
    **kwargs: QueryParam,
) -> Dict[str, Any] | None:
    """Request wrapper that automatically parses JSON response into a python dict.

    This function suppresses NotFoundError and returns None if no results are found.
    """
    try:
        return json.loads(
            get(
                identifier, namespace, domain, operation, "JSON", searchtype, **kwargs
            ).decode()
        )
    except Exception as e:
        logger.info(e)
        return None


def get_compounds(
    identifier: str | int | List[str] | List[int] | List[str | int],
    namespace: str = "cid",
    searchtype: str | None = None,
    **kwargs: QueryParam,
) -> List[str]:
    """Retrieve the specified compound records from PubChem.

    Args:
        identifier: The compound identifier to use as a search query.
        namespace: The identifier type, one of cid, name, smiles, sdf, inchi,
            inchikey or formula.
        searchtype: The advanced search type, one of substructure,
            superstructure or similarity.
        **kwargs: Additional query parameters to pass to the API request.

    Returns:
        List of resolved SMILES strings.
    """
    results = get_json(
        identifier=identifier,
        namespace=namespace,
        domain="compound",
        operation=None,
        searchtype=searchtype,
        **kwargs,
    )

    if results is None:
        raise ValueError("Error retrieving compounds from PubChem")

    if "PC_Compounds" not in results:
        raise ValueError("Error retrieving compounds from PubChem")

    resolved_smiles = []
    for compound in results.get("PC_Compounds", []):
        if compound.get("props", []):
            found_smiles = False
            for prop in compound.get("props", []):
                if prop.get("urn", {}).get("label") == "SMILES":
                    resolved_smiles.append(prop.get("value", {}).get("sval", ""))
                    found_smiles = True
                    break

            if not found_smiles:
                resolved_smiles.append("")
        else:
            resolved_smiles.append("")

    resolved_smiles = ["" if item == DUMMY_SMILES else item for item in resolved_smiles]

    return resolved_smiles


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
        pubchem_smiles = get_compounds(filtered_compound_name_list, "name")
    except Exception as e:
        logger.warning(f"Error with PubChem query: {e}")
        return {}

    pubchem_name_dict = {
        name: smiles if smiles is not None else ""
        for name, smiles in zip(filtered_compound_name_list, pubchem_smiles)
    }

    if len(filtered_compound_name_list) != len(pubchem_name_dict):
        logger.warning(
            f"Mismatching lengths: "
            f"filtered_compound_name_list ({len(filtered_compound_name_list)}), "
            f"pubchem_name_dict ({len(pubchem_name_dict)})"
        )
        return {}

    return pubchem_name_dict
