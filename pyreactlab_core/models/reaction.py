# import libs
from __future__ import annotations

from typing import Any, Dict, Optional, List
from pydantic import BaseModel, Field, computed_field, model_validator
from pythermodb_settings.models import Component
# local imports
from ..core.chem_react import (
    ChemReact,
    ReactionMode,
    PhaseRule
)


class Reaction(BaseModel):
    """
    A class representing a chemical reaction, including its analysis and properties.

    Attributes
    ----------
    name : str
        The name of the reaction.
    reaction : str
        The chemical reaction equation as a string.
    reaction_mode_symbol : Optional[ReactionMode]
        The symbol used to separate reactants and products in a reaction equation.
    analysis : Dict[str, Any]
        A dictionary containing the analysis results of the reaction.

    Properties
    ----------
    symbolic_reaction : str
        The symbolic representation of the balanced reaction.
    symbolic_unbalanced_reaction : str
        The symbolic representation of the unbalanced reaction.
    reactants_names : list[str]
        A list of names of the reactants in the reaction.
    products_names : list[str]
        A list of names of the products in the reaction.
    products : List[Dict[str, Any]]
        A list of dictionaries representing the products of the reaction.
    reactants : List[Dict[str, Any]]
        A list of dictionaries representing the reactants of the reaction.
    reaction_coefficients : Dict[str, float]
        A dictionary of reaction coefficients for each component.
    reaction_stoichiometry : Dict[str, float]
        A dictionary representing the stoichiometry of the reaction.
    reaction_stoichiometry_matrix : list[float]
        A list representing the stoichiometry matrix of the reaction.
    carbon_count : int
        The total number of carbon atoms in the reaction.
    reaction_state : str
        The state of the reaction (e.g., "balanced", "unbalanced").
    reaction_phase : Optional[PhaseRule]
        The phase rule of the reaction, if applicable.
    state_count : Dict[str, int]
        A dictionary counting the states of components in the reaction.
    component_ids : Dict[str, int]
        A dictionary mapping component names to their IDs.
    all_components : list[str]
        A list of all component names involved in the reaction.

    Methods
    -------
    _run_existing_analysis(self) -> Reaction
        Validates and analyzes the reaction after initialization.
    """
    name: str
    reaction: str
    components: Optional[List[Component]] = Field(
        default=None,
        description="A list of Component objects involved in the reaction."
    )
    reaction_mode_symbol: Optional[ReactionMode] = Field(
        default=None,
        description="The symbol used to separate reactants and products in a reaction equation."
    )
    analysis: Dict[str, Any] = Field(
        default_factory=dict,
        description="A dictionary containing the analysis results of the reaction."
    )

    @model_validator(mode="after")
    def _run_existing_analysis(self):
        # NOTE: check reaction mode symbol
        if "<=>" in self.reaction:
            self.reaction_mode_symbol = "<=>"
        elif "=>" in self.reaction:
            self.reaction_mode_symbol = "=>"
        elif "=" in self.reaction:
            self.reaction_mode_symbol = "="
        else:
            raise ValueError(
                f"Invalid reaction format in reaction: {self.reaction}"
            )

        # NOTE: analyze reaction
        util = ChemReact(
            reaction_mode_symbol=self.reaction_mode_symbol,
            components=self.components
        )

        # NOTE: perform analysis
        self.analysis = util.analyze_reaction(
            reaction_pack={
                "name": self.name,
                "reaction": self.reaction
            },
        )
        return self

    @computed_field
    @property
    def symbolic_reaction(self) -> str:
        return self.analysis.get("symbolic_reaction", "")

    @computed_field
    @property
    def symbolic_unbalanced_reaction(self) -> str:
        return self.analysis.get("symbolic_unbalanced_reaction", "")

    @computed_field
    @property
    def reactants_names(self) -> list[str]:
        return self.analysis.get("reactants_names", [])

    @computed_field
    @property
    def products_names(self) -> list[str]:
        return self.analysis.get("products_names", [])

    @computed_field
    @property
    def products(self) -> List[Dict[str, Any]]:
        return self.analysis.get("products", [])

    @computed_field
    @property
    def reactants(self) -> List[Dict[str, Any]]:
        return self.analysis.get("reactants", [])

    @computed_field
    @property
    def reaction_coefficients(self) -> Dict[str, float]:
        return self.analysis.get("reaction_coefficients", {})

    @computed_field
    @property
    def reaction_stoichiometry(self) -> Dict[str, float]:
        return self.analysis.get("reaction_stoichiometry", {})

    @computed_field
    @property
    def reaction_stoichiometry_matrix(self) -> list[float]:
        return self.analysis.get("reaction_stoichiometry_matrix", [])

    @computed_field
    @property
    def carbon_count(self) -> int:
        return self.analysis.get("carbon_count", 0)

    @computed_field
    @property
    def reaction_state(self) -> Dict[str, str]:
        return self.analysis.get("reaction_state", {})

    @computed_field
    @property
    def reaction_phase(self) -> Optional[PhaseRule]:
        return self.analysis.get("reaction_phase", None)

    @computed_field
    @property
    def state_count(self) -> Dict[str, int]:
        return self.analysis.get("state_count", {})

    @computed_field
    @property
    def component_ids(self) -> Dict[str, int]:
        return self.analysis.get("component_ids", {})

    @computed_field
    @property
    def all_components(self) -> list[str]:
        return self.analysis.get("all_components", [])

    @computed_field
    @property
    def available_components(self) -> List[Component]:
        return self.analysis.get("components", [])

    @computed_field
    @property
    def component_checker(self) -> bool:
        return self.analysis.get("component_checker", False)
