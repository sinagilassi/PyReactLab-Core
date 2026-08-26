from __future__ import annotations

from typing import Any, Optional

import numpy as np
from pydantic import BaseModel, Field, computed_field, model_validator
from pythermodb_settings.models import Component
from pythermodb_settings.utils import set_component_id

from ..core.chem_react_balance import ChemReactBalance
from ..core.reaction_network_mapping import (
    component_metadata,
    component_reaction_map,
    reaction_component_map,
)
from ..core.reaction_network_modes import (
    reaction_mode_count,
    reaction_modes,
)
from ..core.reaction_network_phase import (
    components_by_phase,
    multiphase_reactions,
    phases,
    reactions_by_phase,
    single_phase_reactions,
)
from ..core.reaction_network_species import (
    all_species,
    intermediate_species,
    products,
    reactants,
    reaction_ids,
    sink_species,
    source_species,
)
from ..core.reaction_network_stoichiometry import (
    dependent_reaction_ids,
    dependent_reaction_indices,
    independent_reaction_indices,
    independent_reaction_ids,
    independent_stoichiometric_matrix,
    participation_matrix,
    reaction_dependencies,
    stoichiometric_matrix,
    stoichiometric_rank,
)
from .reaction import Reaction


class ReactionNetwork(BaseModel):
    """Structural analysis model for a collection of reactions.

    The stoichiometric matrix uses rows = components/species and columns =
    reactions. Reactants are negative, products are positive, and absent
    species are zero.

    Attributes
    ----------
    name : str
        Name of the reaction network.
    reactions : list[Reaction]
        List of reactions included in the reaction network.
    components : list[Component] | None
        List of chemical components included in the reaction network, the order of which may affect the stoichiometric matrix.

    Notes
    -----
    - The order of components in the `components` list can affect the structure of the stoichiometric matrix.
    - The default component order is determined by component participated in reactions.
    """

    name: str
    reactions: list[Reaction] = Field(
        ...,
        description="Chemical reactions included in the reaction network.",
    )
    components: Optional[list[Component]] = Field(
        default=None,
        description="Chemical components included in the reaction network.",
    )

    @model_validator(mode="after")
    def _validate_network(self):
        if not self.reactions:
            raise ValueError(
                "ReactionNetwork must contain at least one reaction.")

        seen: set[str] = set()
        for reaction in self.reactions:
            if reaction.name in seen:
                raise ValueError(
                    f"Duplicate reaction name detected: {reaction.name}")
            seen.add(reaction.name)

        if self.components is not None:
            if not self.components:
                raise ValueError(
                    "ReactionNetwork components must not be empty when provided."
                )
            self._component_ids_from_components(self.components)

        return self

    @computed_field
    @property
    def reaction_ids(self) -> list[str]:
        return reaction_ids(self.reactions)

    @computed_field
    @property
    def reaction_count(self) -> int:
        return len(self.reactions)

    @computed_field
    @property
    def mapped_components(self) -> list[Any]:
        values: list[Any] = []
        seen: set[str] = set()
        for reaction in self.reactions:
            for component_id, component in reaction.map_components.items():
                if component_id not in seen:
                    seen.add(component_id)
                    values.append(component)
        return values

    @computed_field
    @property
    def component_ids(self) -> list[str]:
        if self.components is None:
            return all_species(self.reactions)
        return self._component_ids_from_components(self.components)

    @computed_field
    @property
    def component_count(self) -> int:
        return len(self.component_ids)

    @computed_field
    @property
    def reactants(self) -> list[str]:
        return reactants(self.reactions)

    @computed_field
    @property
    def products(self) -> list[str]:
        return products(self.reactions)

    @computed_field
    @property
    def all_species(self) -> list[str]:
        return all_species(self.reactions)

    @computed_field
    @property
    def source_species(self) -> list[str]:
        return source_species(self.reactions)

    @computed_field
    @property
    def consumed_species(self) -> list[str]:
        return self.source_species

    @computed_field
    @property
    def sink_species(self) -> list[str]:
        return sink_species(self.reactions)

    @computed_field
    @property
    def produced_species(self) -> list[str]:
        return self.sink_species

    @computed_field
    @property
    def intermediate_species(self) -> list[str]:
        return intermediate_species(self.reactions)

    @computed_field
    @property
    def stoichiometric_matrix(self) -> list[list[float]]:
        return stoichiometric_matrix(
            reactions=self.reactions,
            component_ids=self.component_ids,
        )

    @computed_field
    @property
    def stoichiometric_matrix_dict(self) -> dict[str, list[float]]:
        return {
            component_id: row
            for component_id, row in zip(
                self.component_ids,
                self.stoichiometric_matrix,
            )
        }

    def _component_ids_from_components(
        self,
        components: list[Component],
    ) -> list[str]:
        component_ids = [
            # type: ignore[arg-type]
            set_component_id(component, "Formula-State")
            for component in components
        ]

        duplicate_ids = [
            component_id
            for index, component_id in enumerate(component_ids)
            if component_id in component_ids[:index]
        ]
        if duplicate_ids:
            duplicates = ", ".join(dict.fromkeys(duplicate_ids))
            raise ValueError(f"Duplicate component IDs detected: {duplicates}")

        network_component_ids = set(all_species(self.reactions))
        unknown_ids = [
            component_id
            for component_id in component_ids
            if component_id not in network_component_ids
        ]
        if unknown_ids:
            unknown = ", ".join(unknown_ids)
            raise ValueError(
                f"Components are not part of this network: {unknown}")

        return component_ids

    def stoichiometric_matrix_by_components(
        self,
        components: list[Component] | None = None,
    ) -> list[list[float]]:
        component_ids = (
            self.component_ids
            if components is None
            else self._component_ids_from_components(components)
        )
        return stoichiometric_matrix(
            reactions=self.reactions,
            component_ids=component_ids,
        )

    def stoichiometric_matrix_dict_by_components(
        self,
        components: list[Component] | None = None,
    ) -> dict[str, list[float]]:
        component_ids = (
            self.component_ids
            if components is None
            else self._component_ids_from_components(components)
        )
        matrix = stoichiometric_matrix(
            reactions=self.reactions,
            component_ids=component_ids,
        )
        return {
            component_id: row
            for component_id, row in zip(component_ids, matrix)
        }

    @computed_field
    @property
    def stoichiometric_matrix_shape(self) -> tuple[int, int]:
        return (self.component_count, self.reaction_count)

    @computed_field
    @property
    def stoichiometric_matrix_source(self) -> dict[str, Any]:
        return {
            "components": self.component_ids,
            "reactions": self.reaction_ids,
            "matrix": self.stoichiometric_matrix,
        }

    @computed_field
    @property
    def stoichiometric_rank(self) -> int:
        return stoichiometric_rank(self.stoichiometric_matrix)

    @computed_field
    @property
    def rank(self) -> int:
        return self.stoichiometric_rank

    @computed_field
    @property
    def independent_reactions(self) -> list[str]:
        return independent_reaction_ids(
            matrix=self.stoichiometric_matrix,
            ids=self.reaction_ids,
        )

    @computed_field
    @property
    def independent_reaction_indices(self) -> list[int]:
        return independent_reaction_indices(self.stoichiometric_matrix)

    @computed_field
    @property
    def dependent_reactions(self) -> list[str]:
        return dependent_reaction_ids(
            matrix=self.stoichiometric_matrix,
            ids=self.reaction_ids,
        )

    @computed_field
    @property
    def dependent_reaction_indices(self) -> list[int]:
        return dependent_reaction_indices(
            matrix=self.stoichiometric_matrix,
            independent_indices_value=self.independent_reaction_indices,
        )

    @computed_field
    @property
    def independent_reaction_count(self) -> int:
        return len(self.independent_reactions)

    @computed_field
    @property
    def dependent_reaction_count(self) -> int:
        return len(self.dependent_reactions)

    @computed_field
    @property
    def independent_stoichiometric_matrix(self) -> list[list[float]]:
        return independent_stoichiometric_matrix(
            matrix=self.stoichiometric_matrix,
            independent_indices_value=self.independent_reaction_indices,
        )

    @computed_field
    @property
    def reaction_dependencies(self) -> dict[str, dict[str, float]]:
        return reaction_dependencies(
            matrix=self.stoichiometric_matrix,
            ids=self.reaction_ids,
        )

    @computed_field
    @property
    def reaction_modes(self) -> dict[str, str]:
        return reaction_modes(self.reactions)

    @computed_field
    @property
    def equilibrium_reactions(self) -> list[str]:
        return [
            reaction_id
            for reaction_id, mode in self.reaction_modes.items()
            if mode == "equilibrium"
        ]

    @computed_field
    @property
    def reversible_reactions(self) -> list[str]:
        return [
            reaction_id
            for reaction_id, mode in self.reaction_modes.items()
            if mode == "reversible"
        ]

    @computed_field
    @property
    def irreversible_reactions(self) -> list[str]:
        return [
            reaction_id
            for reaction_id, mode in self.reaction_modes.items()
            if mode == "irreversible"
        ]

    @computed_field
    @property
    def reaction_mode_count(self) -> dict[str, int]:
        return reaction_mode_count(self.reactions)

    @computed_field
    @property
    def component_reaction_map(self) -> dict[str, list[str]]:
        return component_reaction_map(
            reactions=self.reactions,
            component_ids=self.component_ids,
        )

    @computed_field
    @property
    def reaction_component_map(self) -> dict[str, list[str]]:
        return reaction_component_map(
            reactions=self.reactions,
            component_ids=self.component_ids,
        )

    @computed_field
    @property
    def participation_matrix(self) -> list[list[int]]:
        return participation_matrix(self.stoichiometric_matrix)

    @computed_field
    @property
    def phases(self) -> list[str]:
        return phases(self.reactions)

    @computed_field
    @property
    def phase_count(self) -> int:
        return len(self.phases)

    @computed_field
    @property
    def components_by_phase(self) -> dict[str, list[str]]:
        return components_by_phase(
            reactions=self.reactions,
            component_ids=self.component_ids,
        )

    @computed_field
    @property
    def reactions_by_phase(self) -> dict[str, list[str]]:
        return reactions_by_phase(self.reactions)

    @computed_field
    @property
    def single_phase_reactions(self) -> list[str]:
        return single_phase_reactions(self.reactions)

    @computed_field
    @property
    def multiphase_reactions(self) -> list[str]:
        return multiphase_reactions(self.reactions)

    @computed_field
    @property
    def elements(self) -> list[str]:
        parser = ChemReactBalance()
        metadata = component_metadata(self.reactions)
        values: list[str] = []
        for component_id in self.component_ids:
            composition = parser.parse_elemental_composition(
                metadata[component_id]["molecule"]
            )
            values.extend(composition.keys())
        return list(dict.fromkeys(values))

    @computed_field
    @property
    def element_matrix(self) -> list[list[float]]:
        parser = ChemReactBalance()
        metadata = component_metadata(self.reactions)
        rows: list[list[float]] = []
        for element in self.elements:
            row: list[float] = []
            for component_id in self.component_ids:
                composition = parser.parse_elemental_composition(
                    metadata[component_id]["molecule"]
                )
                row.append(float(composition.get(element, 0.0)))
            rows.append(row)
        return rows

    @computed_field
    @property
    def element_balance_matrix(self) -> list[list[float]]:
        if not self.element_matrix:
            return []
        result = np.array(self.element_matrix, dtype=float) @ np.array(
            self.stoichiometric_matrix,
            dtype=float,
        )
        return result.tolist()

    @computed_field
    @property
    def is_element_balanced(self) -> bool:
        return all(reaction.is_element_balanced for reaction in self.reactions)

    @computed_field
    @property
    def unbalanced_element_reactions(self) -> list[str]:
        return [
            reaction.name
            for reaction in self.reactions
            if not reaction.is_element_balanced
        ]

    @computed_field
    @property
    def charge_vector(self) -> list[float]:
        metadata = component_metadata(self.reactions)
        return [
            float(metadata[component_id]["charge"])
            for component_id in self.component_ids
        ]

    @computed_field
    @property
    def charge_balance_vector(self) -> list[float]:
        result = np.array(self.charge_vector, dtype=float) @ np.array(
            self.stoichiometric_matrix,
            dtype=float,
        )
        return result.tolist()

    @computed_field
    @property
    def is_charge_balanced(self) -> bool:
        return all(reaction.is_charge_balanced for reaction in self.reactions)

    @computed_field
    @property
    def unbalanced_charge_reactions(self) -> list[str]:
        return [
            reaction.name
            for reaction in self.reactions
            if not reaction.is_charge_balanced
        ]

    @computed_field
    @property
    def balanced_reactions(self) -> list[str]:
        return [
            reaction.name
            for reaction in self.reactions
            if reaction.is_balanced
        ]

    @computed_field
    @property
    def unbalanced_reactions(self) -> list[str]:
        return [
            reaction.name
            for reaction in self.reactions
            if not reaction.is_balanced
        ]

    @computed_field
    @property
    def is_balanced(self) -> bool:
        return self.is_element_balanced and self.is_charge_balanced

    @computed_field
    @property
    def summary(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "reaction_count": self.reaction_count,
            "component_count": self.component_count,
            "stoichiometric_rank": self.stoichiometric_rank,
            "independent_reaction_count": self.independent_reaction_count,
            "dependent_reaction_count": self.dependent_reaction_count,
            "independent_reactions": self.independent_reactions,
            "dependent_reactions": self.dependent_reactions,
            "source_species": self.source_species,
            "sink_species": self.sink_species,
            "intermediate_species": self.intermediate_species,
            "phases": self.phases,
            "reaction_mode_count": self.reaction_mode_count,
            "is_element_balanced": self.is_element_balanced,
            "is_charge_balanced": self.is_charge_balanced,
            "is_balanced": self.is_balanced,
        }
