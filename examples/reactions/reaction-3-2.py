# import libs
from pyreactlab_core.models.reaction import Reaction
from pythermodb_settings.models import Component
from rich import print

# NOTE: mixed phase reaction
component_caco3 = Component(
    name="Calcium Carbonate",
    formula="CaCO3",
    state="s"
)
component_cao = Component(
    name="Calcium Oxide",
    formula="CaO",
    state="s"
)
component_co2 = Component(
    name="Carbon Dioxide",
    formula="CO2",
    state="g"
)

components_2 = [
    component_caco3,
    component_cao,
    component_co2,
]

reaction_2 = Reaction(
    name="Thermal Decomposition of Calcium Carbonate",
    reaction="CaCO3(s) => CaO(s) + CO2(g)",
    components=components_2,
)

# NOTE: print analysis
print(
    f"[bold underline]Reaction Analysis for: {reaction_2.name}[/bold underline]")
print(f"Reaction: {reaction_2.reaction}")
print(f"Component IDs: {reaction_2.component_ids}")
print(f"Reaction Mode Symbol: {reaction_2.reaction_mode_symbol}")
print(
    f"reaction type: {reaction_2.reaction_type}"
)
print(
    f"Symbolic Unbalanced Reaction: {reaction_2.symbolic_unbalanced_reaction}")
print(f"Symbolic Reaction: {reaction_2.symbolic_reaction}")
print(f"Reactants: {reaction_2.reactants}")
print(f"Products: {reaction_2.products}")
print(f"Reaction Coefficients: {reaction_2.reaction_coefficients}")
print(f"Reaction Stoichiometry: {reaction_2.reaction_stoichiometry}")
print(
    f"Reaction Stoichiometry Matrix: {reaction_2.reaction_stoichiometry_matrix}")
print(
    f"Reaction Stoichiometry Source: {reaction_2.reaction_stoichiometry_source}")
print(f"State Counts: {reaction_2.state_count}")
print(f"Reaction Phase: {reaction_2.reaction_phase}")
print(f"Reaction State: {reaction_2.reaction_state}")
print(f"Carbon Count: {reaction_2.carbon_count}")
print(f"Reactants Names: {reaction_2.reactants_names}")
print(f"Products Names: {reaction_2.products_names}")
print(f"All Components: {reaction_2.all_components}")
print(f"Available components: {reaction_2.available_components}")
print(f"Component Checker: {reaction_2.component_checker}")
print(f"Mapped Components: {reaction_2.map_components}")
