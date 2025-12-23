# import libs
from __future__ import annotations

from typing import Any, Dict, Optional
from pydantic import BaseModel, Field, computed_field, model_validator
# local imports
from ..utils.chem_react import (
    ChemReact,
    ReactionMode,
    PhaseRule
)


class Reaction(BaseModel):
    name: str
    reaction: str
    reaction_mode_symbol: Optional[ReactionMode] = Field(
        default=None,
        description="The symbol used to separate reactants and products in a reaction equation."
    )
    analysis: Dict[str, Any] = Field(default_factory=dict)

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
        util = ChemReact(reaction_mode_symbol=self.reaction_mode_symbol)

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
    def reactants_names(self) -> list[str]:
        return self.analysis.get("reactants_names", [])

    @computed_field
    @property
    def products_names(self) -> list[str]:
        return self.analysis.get("products_names", [])

    @computed_field
    @property
    def products(self) -> list[str]:
        return self.analysis.get("products", [])

    @computed_field
    @property
    def reactants(self) -> list[str]:
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
    def carbon_count(self) -> int:
        return self.analysis.get("carbon_count", 0)

    @computed_field
    @property
    def reaction_state(self) -> str:
        return self.analysis.get("reaction_state", "unknown")

    @computed_field
    @property
    def reaction_phase(self) -> Optional[PhaseRule]:
        return self.analysis.get("reaction_phase", None)

    @computed_field
    @property
    def state_count(self) -> Dict[str, int]:
        return self.analysis.get("state_count", {})
