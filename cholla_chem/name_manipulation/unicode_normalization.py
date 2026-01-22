from typing import List, Tuple, Dict
from cholla_chem.utils.string_utils import clean_strings


def normalize_unicode_and_return_mapping(
    compounds_list: List[str],
) -> Tuple[List[str], Dict[str, str]]:
    """
    Clean a list of strings by replacing certain characters and return a mapping from the original strings to their cleaned versions.

    Args:
        compounds_list (List[str]): A list of strings to clean.

    Returns:
        Tuple[List[str], Dict[str, str]]: A tuple containing the list of cleaned strings and a dictionary mapping the original strings to their cleaned versions.
    """
    cleaned_compounds_list = [clean_strings(ele) for ele in compounds_list]
    cleaned_compounds_dict = {
        k: v for k, v in zip(compounds_list, cleaned_compounds_list)
    }
    return cleaned_compounds_list, cleaned_compounds_dict
