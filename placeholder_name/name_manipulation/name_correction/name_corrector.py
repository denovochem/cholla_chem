from __future__ import annotations
from typing import Dict, List, Optional, Tuple

from placeholder_name.name_manipulation.name_correction.correction_strategies import (
    CorrectionStrategy,
    BracketBalancingStrategy,
    CharacterSubstitutionStrategy,
    LocantCorrectionStrategy,
    PunctuationRestorationStrategy,
)
from placeholder_name.name_manipulation.name_correction.dataclasses import (
    CorrectionCandidate,
    CorrectorConfig,
)
from placeholder_name.name_manipulation.name_correction.scoring import (
    ChemicalNameScorer,
)
from placeholder_name.name_manipulation.name_correction.validators import Validator


class ChemNameCorrector:
    """
    Main class for correcting OCR errors in chemical names.

    This class orchestrates the correction process by:
    1. Applying configured correction strategies
    2. Generating candidate corrections
    3. Scoring candidates
    4. Optionally validating with external tools
    5. Returning ranked results

    Example:
        >>> corrector = ChemNameCorrector()
        >>> results = corrector.correct("2-ch1oropropanoic acid")
        >>> print(results[0].name)
        2-chloropropanoic acid

        >>> # With custom configuration
        >>> config = CorrectorConfig(max_candidates=50)
        >>> corrector = ChemNameCorrector(config)

        >>> # With external validation
        >>> validator = PubChemValidator()
        >>> results = corrector.correct("asprin", validator=validator)

    Attributes:
        config: Configuration for the corrector
        strategies: List of active correction strategies
        scorer: Scoring instance for ranking candidates
    """

    def __init__(
        self,
        config: Optional[CorrectorConfig] = None,
        strategies: Optional[List[CorrectionStrategy]] = None,
    ):
        """
        Initialize the chemical name corrector.

        Args:
            config: Configuration object (uses defaults if None)
            strategies: Custom list of strategies (uses defaults if None)
        """
        self.config = config or CorrectorConfig()
        self.scorer = ChemicalNameScorer(self.config)

        if strategies is not None:
            self.strategies = strategies
        else:
            self.strategies = self._create_default_strategies()

    def _create_default_strategies(self) -> List[CorrectionStrategy]:
        """Create the default set of correction strategies based on config."""
        strategies: List[CorrectionStrategy] = []

        if self.config.enable_locant_correction:
            strategies.append(LocantCorrectionStrategy())

        if self.config.enable_character_substitution:
            char_strategy = CharacterSubstitutionStrategy(
                max_edits=self.config.max_character_substitution_edits
            )
            # Add any custom substitutions
            for orig, replacements in self.config.custom_substitutions.items():
                char_strategy.add_substitution(orig, replacements)
            strategies.append(char_strategy)

        if self.config.enable_punctuation_restoration:
            strategies.append(PunctuationRestorationStrategy())

        if self.config.enable_bracket_balancing:
            strategies.append(BracketBalancingStrategy())

        return strategies

    def add_strategy(self, strategy: CorrectionStrategy) -> None:
        """
        Add a custom correction strategy.

        Args:
            strategy: The strategy to add
        """
        self.strategies.append(strategy)

    def remove_strategy(self, strategy_name: str) -> bool:
        """
        Remove a strategy by name.

        Args:
            strategy_name: Name of the strategy to remove

        Returns:
            True if strategy was found and removed, False otherwise
        """
        for i, strategy in enumerate(self.strategies):
            if strategy.name == strategy_name:
                self.strategies.pop(i)
                return True
        return False

    def correct(
        self,
        name: str,
        validator: Optional[Validator] = None,
        validate_all: bool = False,
    ) -> List[CorrectionCandidate]:
        """
        Correct a chemical name and return ranked candidates.

        Args:
            name: The chemical name to correct
            validator: Optional external validator for candidates
            validate_all: If True, validate all candidates; if False,
                         stop after first valid candidate

        Returns:
            List of CorrectionCandidate objects, sorted by score (descending)
        """
        # Generate all candidates
        candidates = self._generate_all_candidates(name)

        # Remove duplicates while preserving best corrections
        unique_candidates = self._deduplicate_candidates(candidates)

        # Score all candidates
        scored_candidates = [
            self.scorer.score(candidate) for candidate in unique_candidates
        ]

        # Filter by minimum score threshold
        filtered_candidates = [
            c for c in scored_candidates if c.score >= self.config.min_score_threshold
        ]

        # Sort by score (descending)
        sorted_candidates = sorted(
            filtered_candidates, key=lambda c: c.score, reverse=True
        )

        # Limit to max candidates
        limited_candidates = sorted_candidates[: self.config.max_candidates]

        # Optionally validate with external tool
        if validator is not None:
            limited_candidates = self._validate_candidates(
                limited_candidates, validator, validate_all
            )

        return limited_candidates

    def correct_batch(
        self, names: List[str], validator: Optional[Validator] = None
    ) -> Dict[str, List[CorrectionCandidate]]:
        """
        Correct multiple chemical names.

        Args:
            names: List of chemical names to correct
            validator: Optional external validator

        Returns:
            Dictionary mapping original names to their candidates
        """
        results = {}
        for name in names:
            results[name] = self.correct(name, validator)
        return results

    def _generate_all_candidates(self, name: str) -> List[CorrectionCandidate]:
        """Generate candidates from all strategies."""
        candidates: List[CorrectionCandidate] = []

        # Apply each strategy
        current_texts: List[Tuple[str, List[CorrectionCandidate]]] = [
            (name, [])
        ]  # (text, corrections)
        for strategy in self.strategies:
            new_texts = []

            for text, existing_corrections in current_texts:
                new_texts.append((text, existing_corrections))
                for new_text, new_corrections in strategy.generate_candidates(
                    text, self.config
                ):
                    combined_corrections = existing_corrections + new_corrections

                    # Check max corrections limit
                    if (
                        len(combined_corrections)
                        <= self.config.max_corrections_per_candidate
                    ):
                        new_texts.append((new_text, combined_corrections))

                        candidates.append(
                            CorrectionCandidate(
                                name=new_text,
                                original_name=name,
                                corrections=combined_corrections,
                            )
                        )

            # Update current texts for next strategy
            # Limit to prevent explosion
            current_texts = new_texts[: self.config.max_candidates * 2]

        return candidates

    def _deduplicate_candidates(
        self, candidates: List[CorrectionCandidate]
    ) -> List[CorrectionCandidate]:
        """Remove duplicate candidates, keeping the one with fewer corrections."""
        seen: Dict[str, CorrectionCandidate] = {}

        for candidate in candidates:
            if candidate.name not in seen:
                seen[candidate.name] = candidate
            else:
                # Keep the one with fewer corrections
                if candidate.num_corrections < seen[candidate.name].num_corrections:
                    seen[candidate.name] = candidate

        return list(seen.values())

    def _validate_candidates(
        self,
        candidates: List[CorrectionCandidate],
        validator: Validator,
        validate_all: bool,
    ) -> List[CorrectionCandidate]:
        """Validate candidates using external validator."""
        validated_candidates = []

        for candidate in candidates:
            is_valid, result = validator.validate(candidate.name)
            candidate.validated = True

            if is_valid:
                candidate.validation_result = result
                # Boost score for valid candidates
                candidate.score = min(1.0, candidate.score + 0.3)
                validated_candidates.append(candidate)

                if not validate_all:
                    # Add remaining unvalidated candidates
                    remaining_idx = candidates.index(candidate) + 1
                    validated_candidates.extend(candidates[remaining_idx:])
                    break
            else:
                # Lower score for invalid candidates
                candidate.score = max(0.0, candidate.score - 0.2)
                validated_candidates.append(candidate)

        return sorted(validated_candidates, key=lambda c: c.score, reverse=True)

    def get_best_candidate(
        self, name: str, validator: Optional[Validator] = None
    ) -> Optional[CorrectionCandidate]:
        """
        Get the single best correction candidate.

        Args:
            name: Chemical name to correct
            validator: Optional external validator

        Returns:
            Best candidate, or None if no candidates found
        """
        candidates = self.correct(name, validator)
        return candidates[0] if candidates else None

    def explain_corrections(self, candidate: CorrectionCandidate) -> str:
        """
        Generate a human-readable explanation of corrections.

        Args:
            candidate: The candidate to explain

        Returns:
            Multi-line string explaining all corrections
        """
        lines = [
            f"Original: {candidate.original_name}",
            f"Corrected: {candidate.name}",
            f"Score: {candidate.score:.3f}",
            f"Number of corrections: {candidate.num_corrections}",
            "",
            "Score components:",
        ]

        for component, value in candidate.score_components.items():
            lines.append(f"  - {component}: {value:.3f}")

        if candidate.corrections:
            lines.append("")
            lines.append("Corrections applied:")
            for i, correction in enumerate(candidate.corrections, 1):
                lines.append(
                    f"  {i}. [{correction.correction_type.name}] "
                    f"'{correction.original}' → '{correction.replacement}'"
                )
                if correction.description:
                    lines.append(f"     {correction.description}")

        if candidate.validated:
            lines.append("")
            lines.append(f"Validated: {candidate.validation_result or 'No result'}")

        return "\n".join(lines)
