from __future__ import annotations

# import libs
from typing import Any

# locals
from ..configs.constants import ReactionMode
from .reaction_network_species import (
    intermediate_species,
    sink_species,
    source_species,
    unique_preserve_order,
)


# SECTION: legacy compatibility class
class ReactionNetworkAnalysis:
    """
    Backward-compatible utilities for analyzing groups of reactions.
    """

    def analyze_overall_reactions(
        self,
        reactions: list[dict[str, str]],
        reaction_mode_symbol: ReactionMode | None = None,
    ) -> dict[str, list[str]]:
        """
        Analyze raw reaction dictionaries using the Reaction parser.
        """
        # NOTE: local import avoids a module import cycle with ChemReact.
        from ..models.reaction import Reaction

        # SECTION: parse raw reactions through the model layer
        parsed = [
            Reaction(
                name=reaction.get("name", f"R{index + 1}"),
                reaction=reaction["reaction"],
                reaction_mode_symbol=reaction_mode_symbol,
            )
            for index, reaction in enumerate(reactions)
        ]

        # NOTE: keep legacy output keys while using deterministic classifiers.
        return {
            "consumed": source_species(parsed),
            "produced": sink_species(parsed),
            "intermediate": intermediate_species(parsed),
        }

    def analyze_overall_reactions_v2(
        self,
        reactions: dict[str, Any],
    ) -> dict[str, list[str]]:
        """
        Analyze already parsed reaction dictionaries in deterministic order.
        """
        # SECTION: collect parsed reactant and product identifiers
        all_reactants: list[str] = []
        all_products: list[str] = []
        for reaction_value in reactions.values():
            all_reactants.extend(reaction_value["reactants_names"])
            all_products.extend(reaction_value["products_names"])

        # SECTION: classify parsed species
        reactant_values = unique_preserve_order(all_reactants)
        product_values = unique_preserve_order(all_products)
        product_set = set(product_values)
        reactant_set = set(reactant_values)

        return {
            "consumed": [
                species for species in reactant_values if species not in product_set
            ],
            "produced": [
                species for species in product_values if species not in reactant_set
            ],
            "intermediate": [
                species for species in reactant_values if species in product_set
            ],
        }

    def define_component_id(self, reaction_res):
        """
        Define legacy molecule-only component IDs and reaction rows.
        """
        # SECTION: collect unique molecule IDs
        component_list: list[str] = []
        for reaction_data in reaction_res.values():
            component_list.extend(
                item["molecule"] for item in reaction_data["reactants"]
            )
            component_list.extend(
                item["molecule"] for item in reaction_data["products"]
            )

        # NOTE: preserve legacy tuple shape with deterministic ordering.
        component_list = unique_preserve_order(component_list)
        component_dict = {
            component_id: index
            for index, component_id in enumerate(component_list)
        }

        # SECTION: initialize reaction-row dictionaries
        comp_list = [
            {component_id: 0.0 for component_id in component_list}
            for _ in reaction_res
        ]

        # SECTION: fill stoichiometric coefficients
        for row_index, reaction_data in enumerate(reaction_res.values()):
            for reactant in reaction_data["reactants"]:
                comp_list[row_index][reactant["molecule"]] = (
                    -1.0 * float(reactant["coefficient"])
                )
            for product in reaction_data["products"]:
                comp_list[row_index][product["molecule"]] = float(
                    product["coefficient"]
                )

        # SECTION: convert row dictionaries to legacy matrix rows
        comp_coeff = [
            [row[component_id] for component_id in component_list]
            for row in comp_list
        ]

        return component_list, component_dict, comp_list, comp_coeff

    @staticmethod
    def define_component_id_v2(reaction_res):
        """
        Define legacy molecule-state component IDs and reaction rows.
        """
        # SECTION: collect unique molecule-state IDs and state metadata
        component_list: list[str] = []
        component_state_list: list[tuple[str, str, str]] = []
        for reaction_data in reaction_res.values():
            for item in reaction_data["reactants"] + reaction_data["products"]:
                component_list.append(item["molecule_state"])
                component_state_list.append(
                    (item["molecule"], item["state"], item["molecule_state"])
                )

        # NOTE: preserve first-seen order for IDs and metadata.
        component_list = unique_preserve_order(component_list)
        component_state_list = list(dict.fromkeys(component_state_list))
        component_dict = {
            component_id: index
            for index, component_id in enumerate(component_list)
        }

        # SECTION: initialize reaction-row dictionaries
        comp_list = [
            {component_id: 0.0 for component_id in component_list}
            for _ in reaction_res
        ]

        # SECTION: fill stoichiometric coefficients
        for row_index, reaction_data in enumerate(reaction_res.values()):
            for reactant in reaction_data["reactants"]:
                comp_list[row_index][reactant["molecule_state"]] = (
                    -1.0 * float(reactant["coefficient"])
                )
            for product in reaction_data["products"]:
                comp_list[row_index][product["molecule_state"]] = float(
                    product["coefficient"]
                )

        # SECTION: convert row dictionaries to legacy matrix rows
        comp_coeff = [
            [row[component_id] for component_id in component_list]
            for row in comp_list
        ]

        return (
            component_list,
            component_dict,
            comp_list,
            comp_coeff,
            component_state_list,
        )

    def reaction_phase_analysis(
        self,
        reaction_res: dict[str, Any],
    ) -> dict[str, list[str]]:
        """
        Group parsed molecule-state IDs by phase.
        """
        # SECTION: collect component IDs by state
        phase_dict: dict[str, list[str]] = {}
        for reaction_data in reaction_res.values():
            for item in reaction_data["reactants"] + reaction_data["products"]:
                phase_dict.setdefault(item["state"], []).append(
                    item["molecule_state"]
                )

        # NOTE: remove duplicates without changing phase bucket order.
        return {
            phase: unique_preserve_order(component_ids)
            for phase, component_ids in phase_dict.items()
        }
