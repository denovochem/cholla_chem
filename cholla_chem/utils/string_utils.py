import re
from typing import Dict, List

from cholla_chem.utils.constants import (
    ALLOWED_CHARS_WHITELIST,
    COMMON_CHARS_WHITELIST,
    NON_LATIN1_REPLACEMENTS,
)
from cholla_chem.utils.logging_config import logger


def safe_str(x: object) -> str | None:
    """Tries to convert to string, returns none upon exception"""
    try:
        return str(x)
    except Exception:
        return None


def filter_strings_by_whitelist(
    strings: List[str], whitelist: List[str] = COMMON_CHARS_WHITELIST
) -> List[str]:
    """
    Filters a list of strings, returning only those that contain
    only characters from the whitelist.

    :param strings: List of strings to filter.
    :param whitelist: List or set of allowed characters.
    :return: List of strings that contain only whitelisted characters.
    """
    whitelist_set = set(whitelist)  # Convert to set for faster lookup
    return [s for s in strings if set(s).issubset(whitelist_set)]


def clean_strings(
    string: str, chars_to_replace_dict: Dict[str, str] = NON_LATIN1_REPLACEMENTS
) -> str:
    """
    Replaces specified characters in a string.

    Args:
        string: The input string.
        chars_to_replace_dict: A dictionary of characters to replace.

    Returns:
        The cleaned string with the specified characters replaced.
    """
    # sort longest to shortest to avoid partial replacements
    for char in sorted(chars_to_replace_dict, key=len, reverse=True):
        string = string.replace(char, chars_to_replace_dict[char])
    string = string.replace("\n", "")
    string = remove_tags(string)
    string = "".join(char for char in string if char in ALLOWED_CHARS_WHITELIST)
    return string


def is_latin1_compatible(s: str) -> bool:
    """Check if a string is compatible with the latin-1 codec."""
    try:
        s.encode("latin-1")
        return True
    except UnicodeEncodeError:
        logger.info(f"String {s} is not compatible with latin-1 codec")
        return False


def filter_latin1_compatible(strings: List[str]) -> List[str]:
    """Filter a list of strings to only include those compatible with the latin-1 codec."""
    return [s for s in strings if is_latin1_compatible(s)]


def remove_tags(string: str) -> str:
    """Remove HTML tags from a string."""
    pattern = r"<.*?>"
    return re.sub(pattern, "", string)


def cas_check_digit(potential_cas_number: str) -> bool:
    """Check if a string is a valid CAS number."""
    check_digit = potential_cas_number[-1]
    reversed_remaining_digits = potential_cas_number[::-1][1:]
    total = sum(
        (i + 1) * int(digit) for i, digit in enumerate(reversed_remaining_digits)
    )
    if total % 10 == int(check_digit):
        return True
    return False


def is_valid_cas(synonym: str) -> bool:
    """Check if a string is a valid CAS number."""
    synonym = synonym.replace("-", "")
    if not synonym.isdigit():
        return False
    if len(synonym) < 3:
        return False

    return cas_check_digit(synonym)
