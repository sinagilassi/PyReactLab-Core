# import libs
import re
from typing import Dict, List, Optional


# SECTION: ChemReactUtils class
class ChemReactUtils:
    """General-purpose helpers for chemical reaction analysis."""

    # NOTE: supported full phase names
    available_phases = ("gas", "liquid", "aqueous", "solid")

    def __init__(
        self,
        available_phases: tuple[str, ...] | None = None,
    ):
        """
        Initialize general chemical reaction utility settings.
        """
        # SECTION: phase configuration
        # NOTE: child classes can override the supported phase names
        if available_phases is not None:
            self.available_phases = available_phases

    def count_carbon(self, molecule: str, coefficient: float) -> float:
        """
        Count the number of carbon atoms in a molecule.
        """
        try:
            # SECTION: validate inputs
            # NOTE: molecule formula must be text for regex parsing
            if not isinstance(molecule, str):
                raise ValueError("Molecule must be a string.")

            # NOTE: coefficient scales the carbon count
            if not isinstance(coefficient, (int, float)):
                raise ValueError("Coefficient must be an integer or float.")

            # SECTION: carbon symbol matching
            # ! do not count lowercase carbon inside another element symbol
            if re.search(r'C(?![a-z])', molecule):
                # NOTE: multiply atom occurrences by stoichiometric coefficient
                carbon_count = len(re.findall(
                    r'C(?![a-z])', molecule)) * coefficient
                return carbon_count
            else:
                # NOTE: molecule has no carbon atoms
                return 0.0
        except Exception as e:
            raise Exception(
                f"Error counting carbon in molecule '{molecule}': {e}")

    def phase_rule_analysis(self, phase_rule: Optional[str] = None) -> str:
        """
        Analyze the phase rule of a reaction.
        """
        try:
            # SECTION: default phase rule
            # NOTE: empty means component states must be present in the reaction
            if phase_rule is None or phase_rule == 'None':
                return 'empty'

            # SECTION: validate phase rule
            # ? keep this aligned with PhaseRule in chem_react.py
            if phase_rule not in self.available_phases:
                raise ValueError(
                    f"Phase rule must be {', '.join(self.available_phases)}.")

            # SECTION: convert full phase name to reaction state symbol
            if phase_rule == 'gas':
                phase_symbol = 'g'
            elif phase_rule == 'liquid':
                phase_symbol = 'l'
            elif phase_rule == 'aqueous':
                phase_symbol = 'aq'
            elif phase_rule == 'solid':
                phase_symbol = 's'
            else:
                phase_symbol = 'empty'

            # NOTE: return compact state symbol used by parsed components
            return phase_symbol
        except Exception as e:
            raise Exception(f"Error analyzing phase rule: {e}")

    def state_name_set(self, state_set: set) -> List[str]:
        """
        Convert state set to full names.
        """
        try:
            # SECTION: state name mapping
            # NOTE: keys match state symbols parsed from reaction strings
            state_dict = {
                'g': 'gas',
                'l': 'liquid',
                'aq': 'aqueous',
                's': 'solid'
            }

            # NOTE: convert each compact symbol to its full phase name
            return [state_dict[state] for state in state_set]
        except Exception as e:
            raise Exception(f"Error converting state set to full names: {e}")

    def determine_reaction_phase(self, reaction_dict: Dict[str, str]) -> str:
        """
        Determine the phase of a reaction based on component states.
        """
        try:
            # SECTION: collect unique states
            available_states = set(reaction_dict.values())
            # NOTE: convert state symbols before formatting phase text
            state_names = self.state_name_set(available_states)

            # SECTION: determine reaction phase label
            if len(state_names) == 1:
                # NOTE: single-phase reaction
                return f'{state_names[0]}'
            else:
                # NOTE: multi-phase reaction
                return f'{"-".join(state_names)}'
        except Exception as e:
            raise Exception(f"Error determining reaction phase: {e}")

    def count_reaction_states(self, reaction_dict: Dict[str, str]) -> Dict[str, int]:
        """
        Count the number of component states in a reaction.
        """
        try:
            # SECTION: collect component states
            available_states = reaction_dict.values()
            # NOTE: initialize all supported state buckets
            state_count = {
                'g': 0,
                'l': 0,
                'aq': 0,
                's': 0
            }

            # SECTION: count state occurrences
            for state in available_states:
                # ! ignore unsupported states instead of adding new keys
                if state in state_count:
                    state_count[state] += 1

            # NOTE: return counts for every supported state symbol
            return state_count
        except Exception as e:
            raise Exception(f"Error determining reaction phase: {e}")

    # ! reaction types
    def is_reversible(self, reaction_mode_symbol: str) -> bool:
        """
        Check if the reaction is reversible based on the reaction mode symbol.

        Reaction Mode Symbols:
        - `Reversible`: "<=>"
        - `Irreversible`: "=>"
        - `Equilibrium`: "="
        """
        try:
            # SECTION: validate reaction mode symbol
            if reaction_mode_symbol not in ("<=>", "=>", "="):
                raise ValueError(
                    f"Invalid reaction mode symbol: {reaction_mode_symbol}")

            # SECTION: determine reversibility
            return reaction_mode_symbol == "<=>"
        except Exception as e:
            raise Exception(f"Error checking reversibility: {e}")

    def is_irreversible(self, reaction_mode_symbol: str) -> bool:
        """
        Check if the reaction is irreversible based on the reaction mode symbol.

        Reaction Mode Symbols:
        - `Reversible`: "<=>"
        - `Irreversible`: "=>"
        - `Equilibrium`: "="
        """
        try:
            # SECTION: validate reaction mode symbol
            if reaction_mode_symbol not in ("<=>", "=>", "="):
                raise ValueError(
                    f"Invalid reaction mode symbol: {reaction_mode_symbol}")

            # SECTION: determine irreversibility
            return reaction_mode_symbol in ("=>")
        except Exception as e:
            raise Exception(f"Error checking irreversibility: {e}")

    def is_equilibrium(self, reaction_mode_symbol: str) -> bool:
        """
        Check if the reaction is at equilibrium based on the reaction mode symbol.

        Reaction Mode Symbols:
        - `Reversible`: "<=>"
        - `Irreversible`: "=>"
        - `Equilibrium`: "="
        """
        try:
            # SECTION: validate reaction mode symbol
            if reaction_mode_symbol not in ("<=>", "=>", "="):
                raise ValueError(
                    f"Invalid reaction mode symbol: {reaction_mode_symbol}")

            # SECTION: determine equilibrium
            return reaction_mode_symbol == "="
        except Exception as e:
            raise Exception(f"Error checking equilibrium: {e}")

    def reaction_type(self, reaction_mode_symbol: str) -> str:
        """
        Determine the type of reaction based on the reaction mode symbol.

        Reaction Mode Symbols:
        - `Reversible`: "<=>"
        - `Irreversible`: "=>"
        - `Equilibrium`: "="
        """
        try:
            # SECTION: validate reaction mode symbol
            if reaction_mode_symbol not in ("<=>", "=>", "="):
                raise ValueError(
                    f"Invalid reaction mode symbol: {reaction_mode_symbol}")

            # SECTION: determine reaction type
            if self.is_reversible(reaction_mode_symbol):
                return "Reversible"
            elif self.is_irreversible(reaction_mode_symbol):
                return "Irreversible"
            elif self.is_equilibrium(reaction_mode_symbol):
                return "Equilibrium"
            else:
                return "Unknown"
        except Exception as e:
            raise Exception(f"Error determining reaction type: {e}")
