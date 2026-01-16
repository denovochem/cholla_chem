from typing import Optional, Tuple

from typing_extensions import Protocol


class Validator(Protocol):
    """Protocol for external validation of chemical names."""

    def validate(self, name: str) -> Tuple[bool, Optional[str]]:
        """
        Validate a chemical name.

        Args:
            name: Chemical name to validate

        Returns:
            Tuple of (is_valid, result) where result could be SMILES or None
        """
        ...
