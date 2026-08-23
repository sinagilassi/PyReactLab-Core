# import libs
import re
from typing import Any, Dict, List, cast


# SECTION: ReactionNetworkAnalysis class
class ReactionNetworkAnalysis:
    """Utilities that analyze groups of reactions."""

    def analyze_overall_reactions(
            self,
            reactions: List[Dict[str, str]],
            reaction_mode_symbol: str | None = None,
    ) -> Dict[str, List[str]]:
        """
        Analyze a list of chemical reactions and classify species.
        """
        try:
            # SECTION: initialize network species sets
            all_reactants = set()
            all_products = set()

            # SECTION: reaction mode
            # NOTE: prefer explicit method input, fall back to ChemReact state
            reaction_symbol = reaction_mode_symbol
            if reaction_symbol is None:
                reaction_symbol = cast(
                    str,
                    getattr(self, 'reaction_mode_symbol')
                )

            # ! a reaction symbol is required for raw reaction splitting
            if not isinstance(reaction_symbol, str):
                raise ValueError("reaction_mode_symbol must be a string.")

            # SECTION: parse each raw reaction
            for reaction in reactions:
                # NOTE: split equation into reactant and product sides
                sides = reaction['reaction'].split(
                    reaction_symbol.strip()
                )

                # NOTE: simple molecule/state pattern for network classification
                pattern = r'(?:(\d*\.?\d+)\s*)?([A-Z][a-zA-Z0-9]*)\s*(?:\((\w)\))?'

                # SECTION: extract reactant names
                reactants = re.findall(pattern, sides[0])
                reactants = [r[1] for r in reactants]

                # SECTION: extract product names
                products = re.findall(pattern, sides[1])
                products = [p[1] for p in products]

                # NOTE: accumulate species across the whole network
                all_reactants.update(reactants)
                all_products.update(products)

            # SECTION: classify network species
            consumed = list(all_reactants - all_products)
            produced = list(all_products - all_reactants)
            intermediate = list(all_reactants & all_products)

            # NOTE: package classification result
            res = {
                'consumed': consumed,
                'produced': produced,
                'intermediate': intermediate
            }

            # NOTE: return consumed, produced, and intermediate species
            return res
        except Exception as e:
            raise Exception(f"Error analyzing overall reactions: {e}")

    def analyze_overall_reactions_v2(
        self,
        reactions: Dict[str, Any]
    ) -> Dict[str, List[str]]:
        """
        Analyze parsed reactions and classify species.
        """
        try:
            # SECTION: initialize network species sets
            all_reactants = set()
            all_products = set()

            # SECTION: consume already-analyzed reaction dictionaries
            for reaction_name, reaction_value in reactions.items():
                # NOTE: names include molecule-state ids from analyze_reaction
                reactants = reaction_value['reactants_names']
                products = reaction_value['products_names']

                # NOTE: accumulate species across the whole network
                all_reactants.update(reactants)
                all_products.update(products)

            # SECTION: classify network species
            consumed = list(all_reactants - all_products)
            produced = list(all_products - all_reactants)
            intermediate = list(all_reactants & all_products)

            # NOTE: package classification result
            res = {
                'consumed': consumed,
                'produced': produced,
                'intermediate': intermediate
            }

            # NOTE: return consumed, produced, and intermediate species
            return res
        except Exception as e:
            raise Exception(f"Error analyzing overall reactions: {e}")

    def define_component_id(self, reaction_res):
        """
        Define component IDs and stoichiometric coefficients by molecule.
        """
        try:
            # SECTION: collect unique component formulas
            component_list = []

            for item in reaction_res:
                # NOTE: collect reactant formulas
                for reactant in reaction_res[item]['reactants']:
                    component_list.append(reactant['molecule'])
                # NOTE: collect product formulas
                for product in reaction_res[item]['products']:
                    component_list.append(product['molecule'])

            # NOTE: remove duplicate formulas
            component_list = list(set(component_list))

            # SECTION: assign component ids
            component_dict = {}
            for i, item in enumerate(component_list):
                component_dict[item] = i

            # SECTION: initialize stoichiometry rows
            comp_list = [
                {i: 0.0 for i in component_dict.keys()} for _ in range(len(reaction_res))
            ]

            # SECTION: fill stoichiometry rows by reaction
            for j, reaction in enumerate(reaction_res):
                for item in component_dict.keys():
                    # NOTE: reactants have negative coefficients
                    for reactant in reaction_res[reaction]['reactants']:
                        if reactant['molecule'] == item:
                            comp_list[j][item] = -1 * \
                                float(reactant['coefficient'])

                    # NOTE: products have positive coefficients
                    for product in reaction_res[reaction]['products']:
                        if product['molecule'] == item:
                            comp_list[j][item] = float(product['coefficient'])

            # SECTION: convert row dictionaries to matrix rows
            comp_coeff = [
                [comp_list[j][item] for item in component_dict.keys()] for j in range(len(reaction_res))
            ]

            # NOTE: return ids and stoichiometry in legacy tuple format
            return component_list, component_dict, comp_list, comp_coeff
        except Exception as e:
            raise Exception(f"Error defining component ID: {e}")

    @staticmethod
    def define_component_id_v2(reaction_res):
        """
        Define component IDs and stoichiometric coefficients by molecule-state.
        """
        try:
            # SECTION: collect unique molecule-state ids
            component_list = []
            component_state_list = []

            for item in reaction_res:
                # NOTE: collect reactant molecule-state ids
                for reactant in reaction_res[item]['reactants']:
                    component_list.append(reactant['molecule_state'])
                    # NOTE: keep molecule/state metadata alongside ids
                    component_state_list.append(
                        (
                            reactant['molecule'],
                            reactant['state'],
                            reactant['molecule_state']
                        )
                    )
                # NOTE: collect product molecule-state ids
                for product in reaction_res[item]['products']:
                    component_list.append(product['molecule_state'])
                    # NOTE: keep molecule/state metadata alongside ids
                    component_state_list.append(
                        (
                            product['molecule'],
                            product['state'],
                            product['molecule_state']
                        )
                    )

            # NOTE: remove duplicate ids and metadata records
            component_list = list(set(component_list))
            component_state_list = list(set(component_state_list))

            # SECTION: assign component ids
            component_dict = {}

            for i, item in enumerate(component_list):
                component_dict[item] = i

            # SECTION: initialize stoichiometry rows
            comp_list = [
                {i: 0.0 for i in component_dict.keys()} for _ in range(len(reaction_res))
            ]

            # SECTION: fill stoichiometry rows by reaction
            for j, reaction in enumerate(reaction_res):
                for item in component_dict.keys():
                    # NOTE: reactants have negative coefficients
                    for reactant in reaction_res[reaction]['reactants']:
                        if reactant['molecule_state'] == item:
                            comp_list[j][item] = -1 * \
                                float(reactant['coefficient'])

                    # NOTE: products have positive coefficients
                    for product in reaction_res[reaction]['products']:
                        if product['molecule_state'] == item:
                            comp_list[j][item] = float(product['coefficient'])

            # SECTION: convert row dictionaries to matrix rows
            comp_coeff = [
                [comp_list[j][item] for item in component_dict.keys()] for j in range(len(reaction_res))
            ]

            # NOTE: return ids, matrix data, and state metadata in legacy tuple format
            return (
                component_list,
                component_dict,
                comp_list,
                comp_coeff,
                component_state_list
            )
        except Exception as e:
            raise Exception(f"Error defining component ID: {e}")

    def reaction_phase_analysis(
        self,
        reaction_res: Dict[str, Any],
    ):
        """
        Analyze reaction phases and group components by phase.
        """
        try:
            # SECTION: initialize phase buckets
            phase_dict = {
                'g': [],
                'l': [],
                'aq': [],
                's': []
            }

            # SECTION: group participants by phase
            for reaction_name, reaction_data in reaction_res.items():
                # NOTE: group reactants
                for reactant in reaction_data['reactants']:
                    phase = reactant['state']
                    if phase in phase_dict:
                        phase_dict[phase].append(reactant['molecule_state'])
                    else:
                        # ! unsupported state symbol in analyzed reactant
                        raise ValueError(
                            f"Unknown phase '{phase}' for reactant '{reactant['molecule']}'.")

                # NOTE: group products
                for product in reaction_data['products']:
                    phase = product['state']
                    if phase in phase_dict:
                        phase_dict[phase].append(product['molecule_state'])
                    else:
                        # ! unsupported state symbol in analyzed product
                        raise ValueError(
                            f"Unknown phase '{phase}' for product '{product['molecule']}'.")

            # SECTION: remove duplicates from each phase bucket
            for phase in phase_dict:
                phase_dict[phase] = list(set(phase_dict[phase]))

            # NOTE: return molecule-state ids grouped by phase symbol
            return phase_dict
        except Exception as e:
            raise Exception(f"Error analyzing reaction phase: {e}")
