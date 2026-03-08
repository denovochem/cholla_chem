from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field
from enum import Enum, auto
from typing import Any, Dict, List, Optional, Set, Tuple

from rdkit import Chem

from cholla_chem.resolvers.structural_formula_resolver.structural_formula_resolver_tokens import (
    KNOWN_FRAGMENTS,
    ONE_LETTER_ELEMENTS,
    STANDARD_VALENCES,
    TWO_LETTER_ELEMENTS,
    FragmentDefinition,
)
from cholla_chem.utils.logging_config import logger

# =============================================================================
# CONSTANTS AND ELEMENT DATA
# =============================================================================


class BondOrder(Enum):
    """Enumeration of chemical bond orders."""

    SINGLE = 1
    DOUBLE = 2
    TRIPLE = 3
    AROMATIC = 1.5

    def to_smiles_symbol(self) -> str:
        """Convert bond order to SMILES bond symbol."""
        symbols = {
            BondOrder.SINGLE: "",
            BondOrder.DOUBLE: "=",
            BondOrder.TRIPLE: "#",
            BondOrder.AROMATIC: ":",
        }
        return symbols.get(self, "")


def register_fragment(
    pattern: str,
    smiles: str,
    attachment_points: List[int],
    is_aromatic: bool = False,
    description: str = "",
) -> None:
    """
    Register a new fragment pattern for recognition.

    Args:
        pattern: The molecular formula pattern (e.g., 'C6H5')
        smiles: SMILES representation of the fragment
        attachment_points: List of atom indices where bonds can form
        is_aromatic: Whether the fragment is aromatic
        description: Human-readable name

    Example:
        register_fragment('C6H5', 'c1ccccc1', [0], is_aromatic=True, description='phenyl')
    """
    KNOWN_FRAGMENTS[pattern] = FragmentDefinition(
        pattern=pattern,
        smiles=smiles,
        attachment_points=attachment_points,
        is_aromatic=is_aromatic,
        description=description,
    )


# =============================================================================
# DATA STRUCTURES
# =============================================================================


@dataclass
class Atom:
    """Represents an atom in the molecular structure."""

    element: str
    index: int
    implicit_h_count: int = 0
    charge: int = 0

    def __hash__(self) -> int:
        return hash(self.index)

    def __eq__(self, other: object) -> bool:
        return isinstance(other, Atom) and self.index == other.index


@dataclass
class Bond:
    """Represents a bond between two atoms."""

    atom1_idx: int
    atom2_idx: int
    order: BondOrder = BondOrder.SINGLE


class MolecularGraph:
    """
    Represents a molecular structure as a graph.

    Provides methods for building molecules atom-by-atom and bond-by-bond,
    with efficient neighbor lookup for SMILES generation.
    """

    def __init__(self):
        self.atoms: List[Atom] = []
        self.bonds: List[Bond] = []
        self._adjacency: Dict[int, List[Tuple[int, BondOrder]]] = defaultdict(list)

    def add_atom(self, element: str, implicit_h: int = 0, charge: int = 0) -> int:
        """Add an atom and return its index."""
        idx = len(self.atoms)
        self.atoms.append(Atom(element, idx, implicit_h, charge))
        return idx

    def add_bond(
        self, idx1: int, idx2: int, order: BondOrder = BondOrder.SINGLE
    ) -> bool:
        """Add a bond between two atoms. Returns False if bond already exists."""
        if idx1 == idx2:
            return False

        # Check for existing bond
        for neighbor, _ in self._adjacency[idx1]:
            if neighbor == idx2:
                return False

        self.bonds.append(Bond(idx1, idx2, order))
        self._adjacency[idx1].append((idx2, order))
        self._adjacency[idx2].append((idx1, order))
        return True

    def get_neighbors(self, idx: int) -> List[Tuple[int, BondOrder]]:
        """Get all neighbors of an atom with bond orders."""
        return list(self._adjacency[idx])

    def get_atom(self, idx: int) -> Atom:
        """Get atom by index."""
        return self.atoms[idx]

    def atom_count(self) -> int:
        """Return number of atoms."""
        return len(self.atoms)

    def is_empty(self) -> bool:
        """Check if graph has no atoms."""
        return len(self.atoms) == 0

    def copy(self) -> "MolecularGraph":
        """Create a deep copy of this graph."""
        new_graph = MolecularGraph()
        for atom in self.atoms:
            new_graph.add_atom(atom.element, atom.implicit_h_count, atom.charge)
        for bond in self.bonds:
            new_graph.add_bond(bond.atom1_idx, bond.atom2_idx, bond.order)
        return new_graph


# =============================================================================
# TOKENIZER
# =============================================================================


class TokenType(Enum):
    """Types of tokens in structural formula parsing."""

    ELEMENT = auto()
    NUMBER = auto()
    LPAREN = auto()
    RPAREN = auto()
    DOUBLE_BOND = auto()
    TRIPLE_BOND = auto()
    FRAGMENT = auto()
    END = auto()


@dataclass
class Token:
    """A lexical token."""

    type: TokenType
    value: str
    position: int


class Tokenizer:
    """
    Lexical analyzer for structural formulas.

    Converts a formula string into a sequence of tokens for parsing.
    """

    def __init__(self, formula: str):
        self.formula = formula
        self.pos = 0
        self.errors: List[str] = []

    def _try_match_fragment(self) -> Optional[str]:
        """
        Try to match a known fragment pattern at the current position.

        Returns the matched fragment pattern string, or None if no match.
        Fragment patterns are sorted by length (longest first) to ensure
        greedy matching (e.g., C10H7 before C6H5).
        """
        # Sort patterns by length (longest first) for greedy matching
        sorted_patterns = sorted(KNOWN_FRAGMENTS.keys(), key=len, reverse=True)

        for pattern in sorted_patterns:
            if self.formula[self.pos :].startswith(pattern):
                return pattern
        return None

    def tokenize(self) -> List[Token]:
        """Tokenize the formula into a list of tokens."""
        tokens = []

        while self.pos < len(self.formula):
            char = self.formula[self.pos]

            # Skip whitespace
            if char.isspace():
                self.pos += 1
                continue

            # Try to match known fragments first (e.g., C6H5, C10H7)
            # This must come before element matching to catch patterns like C6H5
            fragment = self._try_match_fragment()
            if fragment is not None:
                tokens.append(Token(TokenType.FRAGMENT, fragment, self.pos))
                self.pos += len(fragment)
                continue

            # Try two-letter element first
            if self.pos + 1 < len(self.formula):
                two_char = self.formula[self.pos : self.pos + 2]
                if two_char in TWO_LETTER_ELEMENTS:
                    tokens.append(Token(TokenType.ELEMENT, two_char, self.pos))
                    self.pos += 2
                    continue

            # Single-letter element
            if char in ONE_LETTER_ELEMENTS:
                tokens.append(Token(TokenType.ELEMENT, char, self.pos))
                self.pos += 1
                continue

            # Number
            if char.isdigit():
                start = self.pos
                num_str = ""
                while self.pos < len(self.formula) and self.formula[self.pos].isdigit():
                    num_str += self.formula[self.pos]
                    self.pos += 1
                tokens.append(Token(TokenType.NUMBER, num_str, start))
                continue

            # Parentheses
            if char == "(":
                tokens.append(Token(TokenType.LPAREN, "(", self.pos))
                self.pos += 1
                continue
            if char == ")":
                tokens.append(Token(TokenType.RPAREN, ")", self.pos))
                self.pos += 1
                continue

            # Bond symbols
            if char == "=":
                tokens.append(Token(TokenType.DOUBLE_BOND, "=", self.pos))
                self.pos += 1
                continue
            if char == "#":
                tokens.append(Token(TokenType.TRIPLE_BOND, "#", self.pos))
                self.pos += 1
                continue
            if char == "≡":
                tokens.append(Token(TokenType.TRIPLE_BOND, "≡", self.pos))
                self.pos += 1
                continue

            # Skip hyphens (often used as separators)
            if char == "-":
                self.pos += 1
                continue

            # Unknown character
            self.errors.append(f"Unknown character '{char}' at position {self.pos}")
            self.pos += 1

        tokens.append(Token(TokenType.END, "", self.pos))
        return tokens


# =============================================================================
# PARSER
# =============================================================================


class ParseError(Exception):
    """Raised when parsing fails."""

    pass


@dataclass
class AtomGroup:
    """
    Represents a parsed atom group (e.g., CH3, OH, Cl).

    Attributes:
        central_element: The main element (C, O, N, etc.)
        hydrogen_count: Number of attached hydrogens
        substituents: List of (element, count) for other attached atoms
    """

    central_element: str
    hydrogen_count: int = 0
    substituents: List[Tuple[str, int]] = field(default_factory=list)


@dataclass
class BranchGroup:
    """
    Represents a branch specification like (CH3)2 or (=O).

    Attributes:
        content: The content inside the parentheses. Can be:
                 - AtomGroup for simple single-group branches
                 - int (attachment point index) for complex branches
        multiplier: How many times this branch appears
        bond_order: Bond order for attachment (for things like =O)
        branch_atom_start: Start index of branch atoms (for copying)
        branch_atom_end: End index of branch atoms (exclusive, for copying)
    """

    content: Any
    multiplier: int = 1
    bond_order: BondOrder = BondOrder.SINGLE
    branch_atom_start: Optional[int] = None
    branch_atom_end: Optional[int] = None


class StructuralFormulaParser:
    """
    Parser for condensed structural formulas.

    This parser handles various notational conventions:
    - CH3, CH2, CH, C - carbon with hydrogens
    - (CH3)2, (CH3)3 - multiple identical branches
    - C(O), C(=O) - carbonyl notation
    - Inline branches: CH(CH3)2
    """

    def __init__(self, tokens: List[Token], formula: str = ""):
        self.tokens = tokens
        self.formula = formula
        self.pos = 0
        self.errors: List[str] = []
        self.graph = MolecularGraph()

    def parse(self) -> Optional[MolecularGraph]:
        """Parse tokens into a molecular graph."""
        try:
            self._parse_formula()

            if not self._at_end():
                self.errors.append(f"Unexpected tokens after position {self.pos}")
                return None

            if self.graph.is_empty():
                self.errors.append("Empty molecule")
                return None

            return self.graph

        except ParseError as e:
            self.errors.append(str(e))
            return None
        except Exception as e:
            self.errors.append(f"Parsing error: {e}")
            logger.warning("Unexpected parsing error")
            return None

    # -------------------------------------------------------------------------
    # Token handling utilities
    # -------------------------------------------------------------------------

    def _current(self) -> Token:
        """Get current token."""
        return self.tokens[min(self.pos, len(self.tokens) - 1)]

    def _peek(self, offset: int = 0) -> Token:
        """Look ahead at a token."""
        idx = min(self.pos + offset, len(self.tokens) - 1)
        return self.tokens[idx]

    def _advance(self) -> Token:
        """Consume and return current token."""
        token = self._current()
        if token.type != TokenType.END:
            self.pos += 1
        return token

    def _at_end(self) -> bool:
        """Check if at end of input."""
        return self._current().type == TokenType.END

    def _match(self, *types: TokenType) -> bool:
        """Check if current token matches any given type."""
        return self._current().type in types

    def _consume_number(self) -> int:
        """Consume a number token, return 1 if not present."""
        if self._match(TokenType.NUMBER):
            return int(self._advance().value)
        return 1

    # -------------------------------------------------------------------------
    # Main parsing logic
    # -------------------------------------------------------------------------

    def _parse_formula(self) -> None:
        """
        Parse the entire formula.

        A formula is a sequence of:
        - Prefix branches: (X)n before an atom
        - Atom groups: C, CH, CH2, CH3, OH, etc.
        - Inline branches: (X) after an atom
        - Bond modifiers: =, #
        """
        prev_atom_idx: Optional[int] = None
        pending_branches: List[BranchGroup] = []
        next_bond_order = BondOrder.SINGLE

        while not self._at_end():
            current = self._current()

            # Stop at closing paren (for recursive calls)
            if current.type == TokenType.RPAREN:
                break

            # Bond order modifiers
            if current.type == TokenType.DOUBLE_BOND:
                self._advance()
                next_bond_order = BondOrder.DOUBLE
                continue
            elif current.type == TokenType.TRIPLE_BOND:
                self._advance()
                next_bond_order = BondOrder.TRIPLE
                continue

            # Prefix parenthetical group: (X)n before an atom
            if current.type == TokenType.LPAREN:
                branch = self._parse_parenthetical_branch()
                if branch is not None:
                    pending_branches.append(branch)
                continue

            # Fragment (e.g., C6H5 for phenyl)
            if current.type == TokenType.FRAGMENT:
                fragment_pattern = self._advance().value
                attachment_idx = self._add_fragment_to_graph(fragment_pattern)

                if attachment_idx is None:
                    break

                # Connect to previous atom in chain
                if prev_atom_idx is not None:
                    bond_order_to_use = next_bond_order
                    if bond_order_to_use == BondOrder.SINGLE:
                        bond_order_to_use = self._infer_bond_order(
                            prev_atom_idx, attachment_idx
                        )

                    self.graph.add_bond(
                        prev_atom_idx, attachment_idx, bond_order_to_use
                    )
                    next_bond_order = BondOrder.SINGLE

                # Attach any pending prefix branches to the fragment
                for branch in pending_branches:
                    self._attach_branch(attachment_idx, branch)
                pending_branches.clear()

                prev_atom_idx = attachment_idx

                # Parse any inline branches after the fragment
                while self._match(TokenType.LPAREN):
                    branch = self._parse_parenthetical_branch()
                    if branch is not None:
                        self._attach_branch(attachment_idx, branch)

                continue

            # Atom group
            if current.type == TokenType.ELEMENT:
                atom_idx = self._parse_atom_with_inline_branches(pending_branches)

                if atom_idx is None:
                    break

                # Connect to previous atom in chain
                if prev_atom_idx is not None:
                    # If no explicit bond order was set (i.e., next_bond_order is still SINGLE),
                    # try to infer bond order from valences
                    bond_order_to_use = next_bond_order
                    if bond_order_to_use == BondOrder.SINGLE:
                        bond_order_to_use = self._infer_bond_order(
                            prev_atom_idx, atom_idx
                        )

                    self.graph.add_bond(prev_atom_idx, atom_idx, bond_order_to_use)
                    next_bond_order = BondOrder.SINGLE

                # Attach any pending prefix branches
                for branch in pending_branches:
                    self._attach_branch(atom_idx, branch)
                pending_branches.clear()

                prev_atom_idx = atom_idx
                continue

            # If we get here, we don't know what to do with this token
            self.errors.append(f"Unexpected token: {current}")
            break

        # Any remaining pending branches are an error
        if pending_branches:
            self.errors.append("Prefix branches with no atom to attach to")

    def _infer_bond_order(self, atom1_idx: int, atom2_idx: int) -> BondOrder:
        """
        Infer the most reasonable bond order between two atoms based on their remaining valence.

        Assumes all hydrogens are explicit. Uses STANDARD_VALENCES to compute remaining valence.
        Prefers single > double > triple to avoid overbonding.
        Returns inferred bond order, defaulting to SINGLE if uncertain.
        """

        def get_remaining_valence(atom: Atom) -> int:
            element = atom.element
            if element not in STANDARD_VALENCES:
                return 0  # unknown — can't infer
            possible_valences = STANDARD_VALENCES[element]

            # Choose appropriate valence: prefer common organic valence
            if element == "N":
                # Prefer 3 for neutral N unless already bonded beyond that
                target_valence = 3
                if atom.charge != 0:
                    target_valence = max(possible_valences)
            if element == "S":
                # Prefer 2 for neutral S unless already bonded beyond that
                target_valence = 2
                if atom.charge != 0:
                    target_valence = max(possible_valences)
            elif element in ("C", "P"):
                target_valence = max(possible_valences)
            else:
                target_valence = min(possible_valences)  # O, F, Cl, etc.

            # Sum bond orders (in valence units) from existing bonds
            neighbors = self.graph.get_neighbors(atom.index)
            used_by_bonds = sum(
                self._bond_order_to_valence_units(bo) for _, bo in neighbors
            )

            # Hydrogens count as 1 each
            used_by_h = atom.implicit_h_count

            total_used = used_by_bonds + used_by_h
            remaining = target_valence - total_used
            return max(0, remaining)

        def is_aromatic(atom_idx: int) -> bool:
            for _, bo in self.graph.get_neighbors(atom_idx):
                if bo == BondOrder.AROMATIC:
                    return True
            return False

        if is_aromatic(atom1_idx) or is_aromatic(atom2_idx):
            return BondOrder.SINGLE

        atom1 = self.graph.get_atom(atom1_idx)
        atom2 = self.graph.get_atom(atom2_idx)

        rem1 = get_remaining_valence(atom1)
        rem2 = get_remaining_valence(atom2)

        # Try highest bond order first that both atoms can support
        for order in [BondOrder.TRIPLE, BondOrder.DOUBLE, BondOrder.SINGLE]:
            cost = self._bond_order_to_valence_units(order)
            if rem1 >= cost and rem2 >= cost:
                return order

        return BondOrder.SINGLE  # fallback

    @staticmethod
    def _bond_order_to_valence_units(order: BondOrder) -> int:
        """
        Convert a bond order to the number of valence units it consumes.
        For aromatic bonds, conservatively treat as single (1 unit).
        """
        if order == BondOrder.AROMATIC:
            return 1
        return int(order.value)  # SINGLE=1, DOUBLE=2, TRIPLE=3

    def _parse_parenthetical_branch(self) -> Optional[BranchGroup]:
        """
        Parse a parenthetical expression like (CH3)2 or (=O).

        Returns a BranchGroup describing what was parsed.
        """
        self._advance()  # Consume '('

        # Check for bond order prefix inside parens: (=O), (=S)
        bond_order = BondOrder.SINGLE
        if self._match(TokenType.DOUBLE_BOND):
            self._advance()
            bond_order = BondOrder.DOUBLE
        elif self._match(TokenType.TRIPLE_BOND):
            self._advance()
            bond_order = BondOrder.TRIPLE

        # Parse the content
        content_start_idx = self.graph.atom_count()

        # Simple case: single element like (O), (Cl), (=O)
        if self._match(TokenType.ELEMENT):
            element = self._current().value
            next_tok = self._peek(1)

            # Check if this is a simple single-atom branch
            if next_tok.type == TokenType.RPAREN:
                self._advance()  # consume element
                self._advance()  # consume ')'
                multiplier = self._consume_number()

                # Special case: (O) typically means carbonyl =O
                if element == "O" and bond_order == BondOrder.SINGLE:
                    bond_order = BondOrder.DOUBLE
                elif element == "S" and bond_order == BondOrder.SINGLE:
                    bond_order = BondOrder.DOUBLE

                return BranchGroup(
                    content=AtomGroup(element, 0, []),
                    multiplier=multiplier,
                    bond_order=bond_order,
                )

            # For elements C, N, O, S, try to parse as simple atom group
            # but only accept it if followed by RPAREN
            if element in ("C", "N", "O", "S"):
                saved_pos = self.pos

                try:
                    group = self._parse_simple_atom_group()

                    if self._match(TokenType.RPAREN):
                        self._advance()  # consume ')'
                        multiplier = self._consume_number()
                        return BranchGroup(
                            content=group, multiplier=multiplier, bond_order=bond_order
                        )
                    else:
                        self.pos = saved_pos
                except Exception:
                    self.pos = saved_pos

        # Handle single fragment inside parens: (C6H5), (C2H5)
        if self._match(TokenType.FRAGMENT):
            fragment_pattern = self._advance().value

            if self._match(TokenType.RPAREN):
                self._advance()  # consume ')'
                multiplier = self._consume_number()

                # Build the fragment and record its boundaries
                fragment_start = self.graph.atom_count()
                attachment_idx = self._add_fragment_to_graph(fragment_pattern)
                fragment_end = self.graph.atom_count()

                if attachment_idx is None:
                    return None

                return BranchGroup(
                    content=attachment_idx,
                    multiplier=multiplier,
                    bond_order=bond_order,
                    branch_atom_start=fragment_start,
                    branch_atom_end=fragment_end,
                )
            else:
                # Fragment followed by more content - handle as complex branch
                # First, add the fragment
                fragment_start = self.graph.atom_count()
                inner_prev_atom = self._add_fragment_to_graph(fragment_pattern)
                inner_first_atom = inner_prev_atom
                # Continue parsing below
        else:
            inner_first_atom = None
            inner_prev_atom = None

        inner_bond_order = BondOrder.SINGLE

        while not self._at_end() and not self._match(TokenType.RPAREN):
            if self._match(TokenType.DOUBLE_BOND):
                self._advance()
                inner_bond_order = BondOrder.DOUBLE
                continue
            elif self._match(TokenType.TRIPLE_BOND):
                self._advance()
                inner_bond_order = BondOrder.TRIPLE
                continue

            if self._match(TokenType.ELEMENT):
                group = self._parse_simple_atom_group()
                atom_idx = self._add_atom_group_to_graph(group)

                if inner_first_atom is None:
                    inner_first_atom = atom_idx

                if inner_prev_atom is not None:
                    self.graph.add_bond(inner_prev_atom, atom_idx, inner_bond_order)
                    inner_bond_order = BondOrder.SINGLE

                inner_prev_atom = atom_idx
                continue

            # Handle FRAGMENT tokens inside complex branches
            if self._match(TokenType.FRAGMENT):
                fragment_pattern = self._advance().value
                attachment_idx = self._add_fragment_to_graph(fragment_pattern)

                if attachment_idx is None:
                    break

                if inner_first_atom is None:
                    inner_first_atom = attachment_idx

                if inner_prev_atom is not None:
                    self.graph.add_bond(
                        inner_prev_atom, attachment_idx, inner_bond_order
                    )
                    inner_bond_order = BondOrder.SINGLE

                inner_prev_atom = attachment_idx
                continue

            # Handle nested parentheses
            if self._match(TokenType.LPAREN):
                nested_branch = self._parse_parenthetical_branch()
                if nested_branch is not None and inner_prev_atom is not None:
                    self._attach_branch(inner_prev_atom, nested_branch)
                continue

            break

        if not self._match(TokenType.RPAREN):
            self.errors.append("Unclosed parenthesis")
            return None

        self._advance()  # consume ')'
        multiplier = self._consume_number()

        branch_end = self.graph.atom_count()

        # Use last atom as attachment point (not first!)
        attachment_point = (
            inner_prev_atom if inner_prev_atom is not None else inner_first_atom
        )

        return BranchGroup(
            content=attachment_point,
            multiplier=multiplier,
            bond_order=bond_order,
            branch_atom_start=content_start_idx,
            branch_atom_end=branch_end,
        )

    def _parse_simple_atom_group(self) -> AtomGroup:
        """
        Parse a simple atom group like C, CH, CH2, CH3, OH, NH2, Cl, etc.

        Handles:
        - Element alone: C, O, N, Cl
        - Element with H count: CH3, NH2, OH
        - Element with other substituents: CCl3, CHCl2
        """
        if not self._match(TokenType.ELEMENT):
            raise ParseError("Expected element")

        element = self._advance().value

        # Special case: if element is H, it might be "HC#CH" style
        # In this case, H is attached to what follows
        if element == "H":
            return AtomGroup("H", 0, [])  # Explicit H atom

        hydrogen_count = 0
        substituents: List[Tuple[str, int]] = []

        # Parse attached atoms (H, halogens, etc.)
        while self._match(TokenType.ELEMENT):
            attached = self._current().value

            # Peek ahead: if next is number, it's a count
            if attached == "H":
                self._advance()
                hydrogen_count = self._consume_number()
            elif attached in ("F", "Cl", "Br", "I"):
                # Halogen - can have a count
                self._advance()
                count = self._consume_number()
                substituents.append((attached, count))
            else:
                # Different element - this is the start of the next group
                break

        return AtomGroup(element, hydrogen_count, substituents)

    def _add_atom_group_to_graph(self, group: AtomGroup) -> int:
        """
        Add an AtomGroup to the molecular graph.

        Returns the index of the central atom.
        """
        central_idx = self.graph.add_atom(group.central_element, group.hydrogen_count)

        # Add substituents
        for element, count in group.substituents:
            for _ in range(count):
                # First, add the substituent atom to get its index
                sub_idx = self.graph.add_atom(element, 0)
                # Now infer bond order between central atom and this new atom
                bond_order = self._infer_bond_order(central_idx, sub_idx)
                # Create the bond
                self.graph.add_bond(central_idx, sub_idx, bond_order)

        return central_idx

    def _parse_atom_with_inline_branches(
        self, prefix_branches: List[BranchGroup]
    ) -> Optional[int]:
        """
        Parse an atom group with possible inline branches.

        E.g., "C(O)(CH3)2" or "CH(CH3)"

        Returns the index of the central atom.
        """
        group = self._parse_simple_atom_group()
        central_idx = self._add_atom_group_to_graph(group)

        # Parse inline branches: things like (O), (CH3), (=O) that come after the atom
        while self._match(TokenType.LPAREN):
            branch = self._parse_parenthetical_branch()
            if branch is not None:
                self._attach_branch(central_idx, branch)

        return central_idx

    def _attach_branch(self, attach_to: int, branch: BranchGroup) -> None:
        """
        Attach a branch to an atom.

        Handles multipliers by creating multiple copies of the entire branch structure.
        """
        for i in range(branch.multiplier):
            if isinstance(branch.content, AtomGroup):
                branch_idx = self._add_atom_group_to_graph(branch.content)
                self.graph.add_bond(attach_to, branch_idx, branch.bond_order)
            elif isinstance(branch.content, int):
                if i == 0:
                    # First instance: just bond to existing branch
                    self.graph.add_bond(attach_to, branch.content, branch.bond_order)
                else:
                    # Subsequent instances: copy the entire branch structure
                    if (
                        branch.branch_atom_start is not None
                        and branch.branch_atom_end is not None
                    ):
                        new_attachment = self._copy_branch_structure(
                            branch.branch_atom_start,
                            branch.branch_atom_end,
                            branch.content,  # Original attachment point
                        )
                        if new_attachment is not None:
                            self.graph.add_bond(
                                attach_to, new_attachment, branch.bond_order
                            )
                    else:
                        # Fallback: just copy single atom (for backwards compatibility)
                        if branch.content < len(self.graph.atoms):
                            atom = self.graph.atoms[branch.content]
                            new_idx = self.graph.add_atom(
                                atom.element, atom.implicit_h_count
                            )
                            self.graph.add_bond(attach_to, new_idx, branch.bond_order)

    def _copy_branch_structure(
        self, start_idx: int, end_idx: int, attachment_idx: int
    ) -> Optional[int]:
        """
        Copy a range of atoms and their internal bonds from the graph.

        Args:
            start_idx: Starting atom index (inclusive)
            end_idx: Ending atom index (exclusive)
            attachment_idx: Original attachment point index within the range

        Returns:
            New attachment point index after copying, or None on failure
        """
        if start_idx >= end_idx:
            return None

        if attachment_idx < start_idx or attachment_idx >= end_idx:
            return None

        # Map from old indices to new indices
        idx_map: Dict[int, int] = {}

        # Copy all atoms in the range
        for old_idx in range(start_idx, end_idx):
            atom = self.graph.atoms[old_idx]
            new_idx = self.graph.add_atom(
                atom.element, atom.implicit_h_count, atom.charge
            )
            idx_map[old_idx] = new_idx

        # Copy internal bonds (bonds where both atoms are within the copied range)
        for bond in self.graph.bonds:
            both_in_range = (
                start_idx <= bond.atom1_idx < end_idx
                and start_idx <= bond.atom2_idx < end_idx
            )
            if both_in_range:
                new_atom1 = idx_map[bond.atom1_idx]
                new_atom2 = idx_map[bond.atom2_idx]
                # Only add if this bond doesn't already exist (avoid duplicates)
                self.graph.add_bond(new_atom1, new_atom2, bond.order)

        # Return the new attachment point
        return idx_map.get(attachment_idx)

    def _add_fragment_to_graph(self, fragment_pattern: str) -> Optional[int]:
        """
        Add a known fragment (e.g., C6H5 phenyl) to the molecular graph.
        Builds the fragment directly without requiring RDKit.
        Returns the index of the attachment point atom, or None on failure.
        """
        if fragment_pattern not in KNOWN_FRAGMENTS:
            self.errors.append(f"Unknown fragment pattern: {fragment_pattern}")
            return None

        # Use predefined fragment builders for common patterns
        if fragment_pattern == "C6H5":
            # Phenyl ring: 6 aromatic carbons in a ring
            return self._build_phenyl_ring()
        elif fragment_pattern == "C6H11":
            # Cyclohexyl: 6 sp3 carbons in a ring
            return self._build_cycloalkyl_ring(6)
        elif fragment_pattern == "C5H9":
            # Cyclopentyl: 5 sp3 carbons in a ring
            return self._build_cycloalkyl_ring(5)
        elif fragment_pattern in ["COOH", "HOOC"]:
            # Carboxylic acid
            return self._build_carboxylic_acid()
        elif fragment_pattern == "CONH2":
            # Primary amide
            return self._build_primary_amide()
        elif fragment_pattern == "COOCH3":
            # Methyl ester
            return self._build_ester(n_carbon=1)
        elif fragment_pattern == "COCH3":
            # Methyl ester
            return self._build_aldehyde()
        ## fragment patterns that have multiple attachment points?
        # elif fragment_pattern == "COCH2":
        #     # Methyl ester
        #     return self._build_non_terminal_aldehyde()
        elif fragment_pattern == "COOCH2CH3":
            # Ethyl ester
            return self._build_ester(n_carbon=2)
        elif fragment_pattern == "COCl":
            # Acyl chloride
            return self._build_acid_halide(halide="Cl")
        elif fragment_pattern == "COBr":
            # Acyl bromide
            return self._build_acid_halide(halide="Br")
        elif fragment_pattern == "COI":
            # Acyl iodide
            return self._build_acid_halide(halide="I")
        elif fragment_pattern == "C2H5":
            return self._build_alkyl_chain(2)
        elif fragment_pattern == "C3H7":
            return self._build_alkyl_chain(3)
        elif fragment_pattern == "C4H9":
            return self._build_alkyl_chain(4)
        elif fragment_pattern == "NO2":
            return self._build_nitro()
        elif fragment_pattern == "O2N":
            return self._build_nitro()

        return None

    def _build_phenyl_ring(self) -> int:
        """Build a phenyl (benzene) ring and return the attachment point index."""
        # Add 6 aromatic carbons
        indices = []
        for _ in range(6):
            idx = self.graph.add_atom("C", 0)
            indices.append(idx)

        # Connect them in a ring with aromatic bonds
        for i in range(6):
            next_i = (i + 1) % 6
            self.graph.add_bond(indices[i], indices[next_i], BondOrder.AROMATIC)

        # Return the first carbon as the attachment point
        return indices[0]

    def _build_cycloalkyl_ring(self, size: int) -> int:
        """Build a cycloalkyl ring of given size and return the attachment point."""
        indices = []
        for _ in range(size):
            idx = self.graph.add_atom("C", 0)
            indices.append(idx)

        # Connect them in a ring with single bonds
        for i in range(size):
            next_i = (i + 1) % size
            self.graph.add_bond(indices[i], indices[next_i], BondOrder.SINGLE)

        return indices[0]

    def _build_alkyl_chain(self, length: int) -> int:
        """
        Build a straight-chain alkyl group and return the attachment point index.

        The chain is built as CH3-(CH2)n- where the last CH2 is the attachment point.
        """
        if length < 1:
            return -1

        indices = []
        for i in range(length):
            # First carbon (terminal): 3 hydrogens
            # Other carbons: 2 hydrogens each
            h_count = 3 if i == 0 else 2
            idx = self.graph.add_atom("C", h_count)
            indices.append(idx)

        # Bond the carbons in sequence
        for i in range(length - 1):
            self.graph.add_bond(indices[i], indices[i + 1], BondOrder.SINGLE)

        # Attachment point is the last carbon
        return indices[-1]

    def _build_carboxylic_acid(self) -> int:
        """Build a carboxylic acid functional group"""
        idx1 = self.graph.add_atom("C", 0)
        idx2 = self.graph.add_atom("O", 0)
        idx3 = self.graph.add_atom("O", 1)

        self.graph.add_bond(idx1, idx2, BondOrder.DOUBLE)
        self.graph.add_bond(idx1, idx3, BondOrder.SINGLE)

        return idx1

    def _build_nitro(self) -> int:
        """Build a carboxylic acid functional group"""
        idx1 = self.graph.add_atom("N", 1)
        idx2 = self.graph.add_atom("O", 0)
        idx3 = self.graph.add_atom("O", -1)

        self.graph.add_bond(idx1, idx2, BondOrder.DOUBLE)
        self.graph.add_bond(idx1, idx3, BondOrder.DOUBLE)

        return idx1

    def _build_primary_amide(self) -> int:
        """Build a carboxylic acid functional group"""
        idx1 = self.graph.add_atom("C", 0)
        idx2 = self.graph.add_atom("O", 0)
        idx3 = self.graph.add_atom("N", 0)

        self.graph.add_bond(idx1, idx2, BondOrder.DOUBLE)
        self.graph.add_bond(idx1, idx3, BondOrder.SINGLE)

        return idx1

    def _build_ester(self, n_carbon: int = 1) -> int:
        """Build a carboxylic acid functional group"""
        idx1 = self.graph.add_atom("C", 0)
        idx2 = self.graph.add_atom("O", 0)
        idx3 = self.graph.add_atom("O", 0)

        self.graph.add_bond(idx1, idx2, BondOrder.DOUBLE)
        self.graph.add_bond(idx1, idx3, BondOrder.SINGLE)

        prev_idx = idx3
        for _ in range(n_carbon):
            next_idx = self.graph.add_atom("C", 2)
            self.graph.add_bond(prev_idx, next_idx, BondOrder.SINGLE)
            prev_idx = next_idx

        return idx1

    def _build_aldehyde(self) -> int:
        """Build a carboxylic acid functional group"""
        idx1 = self.graph.add_atom("C", 0)
        idx2 = self.graph.add_atom("O", 0)
        idx3 = self.graph.add_atom("C", 0)

        self.graph.add_bond(idx1, idx2, BondOrder.DOUBLE)
        self.graph.add_bond(idx1, idx3, BondOrder.SINGLE)

        return idx1

    def _build_non_terminal_aldehyde(self) -> List[int]:
        """Build a carboxylic acid functional group"""
        idx1 = self.graph.add_atom("C", 0)
        idx2 = self.graph.add_atom("O", 0)
        idx3 = self.graph.add_atom("C", 0)

        self.graph.add_bond(idx1, idx2, BondOrder.DOUBLE)
        self.graph.add_bond(idx1, idx3, BondOrder.SINGLE)

        return [idx1, idx3]

    def _build_acid_halide(self, halide: str = "Cl") -> int:
        """Build an ester functional group"""
        idx1 = self.graph.add_atom("C", 0)
        idx2 = self.graph.add_atom("O", 0)
        idx3 = self.graph.add_atom(halide, 0)

        self.graph.add_bond(idx1, idx2, BondOrder.DOUBLE)
        self.graph.add_bond(idx1, idx3, BondOrder.SINGLE)

        return idx1


# =============================================================================
# SMILES GENERATOR
# =============================================================================


class SMILESGenerator:
    """
    Generates SMILES strings from molecular graphs using RDKit.

    This is a drop-in replacement for the custom DFS-based SMILES generator.
    It builds an RDKit molecule from the molecular graph and uses RDKit's
    SMILES generation for more robust and standardized output.
    """

    def __init__(self, graph: MolecularGraph):
        """
        Initialize the generator with a molecular graph.

        Args:
            graph: The molecular graph to convert to SMILES.
        """
        self.graph = graph

    def generate(self) -> str:
        """
        Generate a SMILES string from the molecular graph.

        Returns:
            Canonical SMILES string, or empty string if generation fails.
        """
        if self.graph.is_empty():
            return ""

        try:
            mol = self._build_rdkit_mol()
            if mol is None:
                return ""

            # Attempt sanitization with fallback strategies
            if not self._sanitize_mol(mol):
                logger.warning("Failed to sanitize molecule during SMILES generation")
                return ""

            # Generate canonical SMILES
            return Chem.MolToSmiles(mol, canonical=True)

        except Exception as e:
            logger.warning(f"SMILES generation failed: {e}")
            return ""

    def _build_rdkit_mol(self) -> Optional[Chem.RWMol]:
        """
        Build an RDKit editable molecule (RWMol) from the molecular graph.

        Atoms are added in order to preserve index correspondence with the
        original graph. Aromatic atoms are marked based on their participation
        in aromatic bonds.

        Returns:
            RWMol object with atoms and bonds added, or None on failure.
        """
        mol = Chem.RWMol()

        # Identify aromatic atoms based on their bonds
        aromatic_atom_indices = self._find_aromatic_atoms()

        # Add atoms in order (preserving graph indices)
        for atom in self.graph.atoms:
            rd_atom = Chem.Atom(atom.element)
            rd_atom.SetFormalCharge(atom.charge)

            # Mark atoms in aromatic systems
            if atom.index in aromatic_atom_indices:
                rd_atom.SetIsAromatic(True)

            # Note: We don't explicitly set hydrogen counts here.
            # RDKit will calculate implicit hydrogens based on valence
            # after sanitization, which matches the original behavior.

            mol.AddAtom(rd_atom)

        # Add bonds with appropriate bond types
        for bond in self.graph.bonds:
            bond_type = self._convert_bond_order(bond.order)
            mol.AddBond(bond.atom1_idx, bond.atom2_idx, bond_type)

        return mol

    def _sanitize_mol(self, mol: Chem.RWMol) -> bool:
        """
        Attempt to sanitize the molecule using various strategies.

        Tries full sanitization first, then falls back to partial sanitization
        if that fails (e.g., due to aromaticity perception conflicts with
        explicitly set aromatic flags).

        Args:
            mol: RDKit molecule to sanitize.

        Returns:
            True if sanitization succeeded, False otherwise.
        """
        # Strategy 1: Full sanitization
        try:
            Chem.SanitizeMol(mol)
            return True
        except Exception:
            pass

        # Strategy 2: Skip aromaticity perception since we set it explicitly
        try:
            Chem.SanitizeMol(
                mol,
                sanitizeOps=(
                    Chem.SanitizeFlags.SANITIZE_ALL
                    ^ Chem.SanitizeFlags.SANITIZE_SETAROMATICITY
                ),
            )
            return True
        except Exception:
            pass

        # Strategy 3: Minimal sanitization as last resort
        try:
            Chem.SanitizeMol(
                mol,
                sanitizeOps=(
                    Chem.SanitizeFlags.SANITIZE_FINDRADICALS
                    | Chem.SanitizeFlags.SANITIZE_SETCONJUGATION
                    | Chem.SanitizeFlags.SANITIZE_SETHYBRIDIZATION
                ),
            )
            return True
        except Exception:
            pass

        return False

    def _find_aromatic_atoms(self) -> Set[int]:
        """
        Identify atoms that participate in aromatic bonds.

        Returns:
            Set of atom indices that have at least one aromatic bond.
        """
        aromatic = set()
        for bond in self.graph.bonds:
            if bond.order == BondOrder.AROMATIC:
                aromatic.add(bond.atom1_idx)
                aromatic.add(bond.atom2_idx)
        return aromatic

    @staticmethod
    def _convert_bond_order(order: BondOrder) -> Chem.BondType:
        """
        Convert internal BondOrder enum to RDKit BondType.

        Args:
            order: The internal bond order representation.

        Returns:
            Corresponding RDKit BondType.
        """
        conversion_map = {
            BondOrder.SINGLE: Chem.BondType.SINGLE,
            BondOrder.DOUBLE: Chem.BondType.DOUBLE,
            BondOrder.TRIPLE: Chem.BondType.TRIPLE,
            BondOrder.AROMATIC: Chem.BondType.AROMATIC,
        }
        return conversion_map.get(order, Chem.BondType.SINGLE)


# =============================================================================
# MAIN CONVERTER
# =============================================================================


class StructuralFormulaConverter:
    """
    Main interface for converting structural formulas to SMILES.

    This converter combines:
    1. Pattern matching for known compounds (most reliable)
    2. Parsing-based conversion for general cases
    3. Validation to ensure results are reasonable

    Usage:
        converter = StructuralFormulaConverter()
        smiles = converter.convert("CH3CH2OH")  # Returns "CCO"

        # Check for errors
        if smiles is None:
            print(converter.get_errors())
    """

    def __init__(self, strict_mode: bool = True, use_registry: bool = True):
        """
        Initialize the converter.

        Args:
            strict_mode: If True, return None for ambiguous/uncertain results
            use_registry: If True, check pattern registry first
        """
        self.strict_mode = strict_mode
        self.use_registry = use_registry
        self.last_errors: List[str] = []

    def convert(self, formula: str) -> str:
        """
        Convert a structural formula to SMILES.

        Args:
            formula: Condensed structural formula

        Returns:
            SMILES string if successful, empty string otherwise
        """
        self.last_errors = []

        # Preprocess
        if not formula:
            self.last_errors.append("Empty formula")
            return ""

        # Tokenize
        tokenizer = Tokenizer(formula)
        tokens = tokenizer.tokenize()

        if tokenizer.errors:
            self.last_errors.extend(tokenizer.errors)
            if self.strict_mode:
                return ""

        # Parse
        parser = StructuralFormulaParser(tokens, formula)
        graph = parser.parse()

        if parser.errors:
            self.last_errors.extend(parser.errors)

        if graph is None:
            return ""

        # Generate SMILES
        generator = SMILESGenerator(graph)
        smiles = generator.generate()

        if not smiles:
            self.last_errors.append("Failed to generate SMILES")
            return ""

        if not self._validate_smiles(smiles):
            self.last_errors.append(
                f"Generated SMILES '{smiles}' failed RDKit validation"
            )
            return ""

        return smiles

    def _validate_smiles(self, smiles: str) -> bool:
        """Validate SMILES using RDKit."""
        try:
            from rdkit import Chem

            mol = Chem.MolFromSmiles(smiles)
            return mol is not None
        except Exception:
            return False

    def get_errors(self) -> List[str]:
        """Get errors from last conversion."""
        return list(self.last_errors)

    def batch_convert(self, formulas: List[str]) -> Dict[str, str]:
        """Convert multiple formulas."""
        return {f: self.convert(f) for f in formulas}


def name_to_smiles_structural_formula(
    formulas: List[str], strict: bool = True
) -> Dict[str, str]:
    """Convert multiple formulas to SMILES."""
    converter = StructuralFormulaConverter(strict_mode=strict)
    return converter.batch_convert(formulas)
