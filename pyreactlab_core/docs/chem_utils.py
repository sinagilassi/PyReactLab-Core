# import libs
import logging
from typing import List, Optional
# locals
from ..models.reaction import Reaction

# NOTE: configure logger
logger = logging.getLogger(__name__)


def build_stoichiometry_matrix(reactions: List[Reaction]):
    '''
    Build stoichiometry matrix for reactions

    Parameters
    ----------
    reactions : List[Reaction]
        List of Reaction instances

    Returns
    -------
    component_list: list
        component list
    component_dict: dict
        component dict
    comp_list: list
        component list
    comp_coeff: list
        component coefficient
    component_state_list: list
        component state list
    '''
    try:
        # SECTION: extract reaction results
        # NOTE: reaction num
        reaction_num = len(reactions)

        # NOTE: component list
        component_list = []

        # NOTE: component state list
        component_state_list = []

        # SECTION: Iterate over reactions and extract reactants and products
        for item in reactions:
            # get components
            _components = item.all_components
            # store
            component_list.extend(_components)

        # remove duplicate
        component_list = list(set(component_list))

        # component id: key, value
        component_dict = {}

        # loop over component list
        for i, item in enumerate(component_list):
            component_dict[item] = i

        # SECTION: Initialize the component list
        comp_list = [
            {i: 0.0 for i in component_dict.keys()} for _ in range(reaction_num)
        ]

        # SECTION: Iterate over reactions and components
        for j, reaction in enumerate(reactions):
            for item in component_dict.keys():
                # NOTE: Check reactants
                for reactant in reactions[j].reactants:
                    # matching state
                    if reactant['molecule_state'] == item:
                        comp_list[j][item] = -1 * \
                            float(reactant['coefficient'])

                    # >> component state list
                    component_state_list.append(
                        (
                            reactant['molecule'],
                            reactant['state'],
                            reactant['molecule_state']
                        )
                    )

                # NOTE: Check products
                for product in reactions[j].products:
                    # matching state
                    if product['molecule_state'] == item:
                        comp_list[j][item] = float(product['coefficient'])

                    # >> component state list
                    component_state_list.append(
                        (
                            product['molecule'],
                            product['state'],
                            product['molecule_state']
                        )
                    )

        # Convert comp_list to comp_matrix
        comp_coeff = [
            [comp_list[j][item] for item in component_dict.keys()] for j in range(reaction_num)
        ]

        # res
        return {
            "component_list": component_list,
            "component_dict": component_dict,
            "component_state_list": component_state_list,
            "comp_list": comp_list,
            "comp_coeff": comp_coeff,
        }
    except Exception as e:
        raise Exception(f"Error defining component ID: {e}")
