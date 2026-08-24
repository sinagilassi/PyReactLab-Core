# import libs
import logging
from typing import Any
from rich import print

# NOTE: setup logger
logger = logging.getLogger(__name__)

# ::: print reaction analysis


def print_reaction_analysis(reaction: Any) -> None:
    """Print a complete formatted analysis for a reaction.

    Parameters
    ----------
    reaction : Any
        A reaction object exposing the analysis attributes used by this
        function, including reaction details, stoichiometry, component data,
        carbon counts, and charge totals.

    Returns
    -------
    None
        Prints each analysis label and its value using Rich formatting.

    Notes
    -----
    The function assumes that the `reaction` object has the following attributes:
    - name
    - reaction
    - component_ids
    - reaction_mode_symbol
    - reaction_type
    - symbolic_unbalanced_reaction
    - symbolic_reaction
    - reactants
    - products
    - reaction_coefficients
    - reaction_stoichiometry
    - reaction_stoichiometry_matrix
    - reaction_stoichiometry_source
    - state_count
    - reaction_phase
    - reaction_state
    - carbon_count
    - net_carbon_count
    - reactants_names
    - products_names
    - all_components
    - available_components
    - component_checker
    - map_components
    - charge_count
    - total_reactant_charge
    - total_product_charge
    - net_charge

    """
    print(
        f"[bold underline]Reaction Analysis for: {reaction.name}[/bold underline]")
    print("Reaction: ")
    print(reaction.reaction)
    print("Component IDs: ")
    print(reaction.component_ids)
    print("Reaction Mode Symbol: ")
    print(reaction.reaction_mode_symbol)
    print("Reaction Type: ")
    print(reaction.reaction_type)
    print("Symbolic Unbalanced Reaction: ")
    print(reaction.symbolic_unbalanced_reaction)
    print("Symbolic Reaction: ")
    print(reaction.symbolic_reaction)
    print("Reactants: ")
    print(reaction.reactants)
    print("Products: ")
    print(reaction.products)
    print("Reaction Coefficients: ")
    print(reaction.reaction_coefficients)
    print("Reaction Stoichiometry: ")
    print(reaction.reaction_stoichiometry)
    print("Reaction Stoichiometry Matrix: ")
    print(reaction.reaction_stoichiometry_matrix)
    print("Reaction Stoichiometry Source: ")
    print(reaction.reaction_stoichiometry_source)
    print("State Counts: ")
    print(reaction.state_count)
    print("Reaction Phase: ")
    print(reaction.reaction_phase)
    print("Reaction State: ")
    print(reaction.reaction_state)
    print("Carbon Count: ")
    print(reaction.carbon_count)
    print("Carbon counts: ")
    print(reaction.net_carbon_count)
    print("Reactants Names: ")
    print(reaction.reactants_names)
    print("Products Names: ")
    print(reaction.products_names)
    print("All Components: ")
    print(reaction.all_components)
    print("Available components: ")
    print(reaction.available_components)
    print("Component Checker: ")
    print(reaction.component_checker)
    print("Mapped Components: ")
    print(reaction.map_components)
    print("Charge Count: ")
    print(reaction.charge_count)
    print("total reactant charge: ")
    print(reaction.total_reactant_charge)
    print("total product charge: ")
    print(reaction.total_product_charge)
    print("net charge: ")
    print(reaction.net_charge)
