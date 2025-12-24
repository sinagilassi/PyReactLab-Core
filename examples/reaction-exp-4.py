# import libs
from pyreactlab_core import rxn, rxn_stoichiometry, rxns_stoichiometry
from rich import print


# NOTE: define reaction string
reaction_1 = "CO2(g) + 3H2(g) => CH3OH(g) + H2O(g)"
name_1 = "CO2 Hydrogenation to Methanol"

# second reaction
reaction_2 = "C2H4(g) + H2(g) => C2H6(g)"
name_2 = "Ethylene Hydrogenation to Ethane"


# NOTE: create reaction instance
rxn_1 = rxn(
    reaction_str=reaction_1,
    name=name_1
)
# NOTE: display analysis results
if rxn_1 is None:
    raise ValueError("Failed to create Reaction instance.")

rxn_2 = rxn(
    reaction_str=reaction_2,
    name=name_2
)
# NOTE: display analysis results
if rxn_2 is None:
    raise ValueError("Failed to create Reaction instance.")

# SECTION: display reaction analysis results
print(f"[bold green]Reaction Name:[/bold green] {rxn_1.name}")
print(
    f"[bold green]Symbolic Reaction:[/bold green] {rxn_1.symbolic_reaction}")
print(f"[bold green]Reactants:[/bold green] {rxn_1.reactants}")
print(f"[bold green]Products:[/bold green] {rxn_1.products}")
print(
    f"[bold green]Reaction Coefficients:[/bold green] {rxn_1.reaction_coefficients}")
print(
    f"[bold green]Reaction Stoichiometry:[/bold green] {rxn_1.reaction_stoichiometry}")
print(
    f"[bold green]Reaction Stoichiometry Matrix:[/bold green] {rxn_1.reaction_stoichiometry_matrix}")
print(
    f"[bold green]Carbon Count:[/bold green] {rxn_1.carbon_count}")
print(
    f"[bold green]Reaction State:[/bold green] {rxn_1.reaction_state}")


# SECTION: Get reaction stoichiometry matrix
stoichiometry_result = rxn_stoichiometry(
    reaction=rxn_1,
)
if stoichiometry_result is None:
    raise ValueError("Failed to retrieve reaction stoichiometry matrix.")
print(
    f"[bold blue]Stoichiometry Matrix from rxn_stoichiometry():[/bold blue] ")
print(stoichiometry_result)

# SECTION: Get stoichiometry matrices for multiple reactions
reactions_list = [rxn_1, rxn_2]
stoichiometry_matrices = rxns_stoichiometry(
    reactions=reactions_list,
)
if stoichiometry_matrices is None:
    raise ValueError(
        "Failed to retrieve stoichiometry matrices for reactions.")

print(stoichiometry_matrices)
