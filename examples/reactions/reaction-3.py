# import libs
from pyreactlab_core.models.reaction import Reaction
from pythermodb_settings.models import Component
from rich import print

# NOTE: gas phase reaction
component_h2 = Component(
    name="Hydrogen",
    formula="H2",
    state="g"
)
component_o2 = Component(
    name="Oxygen",
    formula="O2",
    state="g"
)
component_h2o_g = Component(
    name="Water",
    formula="H2O",
    state="g"
)

components_1 = [
    component_h2,
    component_o2,
    component_h2o_g,
]

reaction_1 = Reaction(
    name="Combustion of Hydrogen",
    reaction="2H2(g) + O2(g) => 2H2O(g)",
    components=components_1,
)

# NOTE: print analysis
print(
    f"[bold underline]Reaction Analysis for: {reaction_1.name}[/bold underline]")
print(f"Reaction: {reaction_1.reaction}")
print(f"Component IDs: {reaction_1.component_ids}")
print(f"Reaction Mode Symbol: {reaction_1.reaction_mode_symbol}")
print(
    f"reaction type: {reaction_1.reaction_type}"
)
print(
    f"Symbolic Unbalanced Reaction: {reaction_1.symbolic_unbalanced_reaction}")
print(f"Symbolic Reaction: {reaction_1.symbolic_reaction}")
print(f"Reactants: {reaction_1.reactants}")
print(f"Products: {reaction_1.products}")
print(f"Reaction Coefficients: {reaction_1.reaction_coefficients}")
print(f"Reaction Stoichiometry: {reaction_1.reaction_stoichiometry}")
print(
    f"Reaction Stoichiometry Matrix: {reaction_1.reaction_stoichiometry_matrix}")
print(
    f"Reaction Stoichiometry Source: {reaction_1.reaction_stoichiometry_source}")
print(f"State Counts: {reaction_1.state_count}")
print(f"Reaction Phase: {reaction_1.reaction_phase}")
print(f"Reaction State: {reaction_1.reaction_state}")
print(f"Carbon Count: {reaction_1.carbon_count}")
print(f"Reactants Names: {reaction_1.reactants_names}")
print(f"Products Names: {reaction_1.products_names}")
print(f"All Components: {reaction_1.all_components}")
print(f"Available components: {reaction_1.available_components}")
print(f"Component Checker: {reaction_1.component_checker}")
print(f"Mapped Components: {reaction_1.map_components}")
