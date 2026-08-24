# import libs
from pyreactlab_core.models.reaction import Reaction
from pythermodb_settings.models import Component
from rich import print

# NOTE: aqueous neutralization reaction
component_h_ion = Component(
    name="Hydrogen Ion",
    formula="H",
    state="aq"
)
component_oh_ion = Component(
    name="Hydroxide Ion",
    formula="OH",
    state="aq"
)
component_h2o_l = Component(
    name="Water",
    formula="H2O",
    state="l"
)

components_3 = [
    component_h_ion,
    component_oh_ion,
    component_h2o_l,
]

reaction_3 = Reaction(
    name="Acid-Base Neutralization",
    reaction="H{+}(aq) + OH{-}(aq) => H2O(l)",
    components=components_3,
)

# NOTE: print analysis
print(
    f"[bold underline]Reaction Analysis for: {reaction_3.name}[/bold underline]")
print(f"Reaction: {reaction_3.reaction}")
print(f"Component IDs: {reaction_3.component_ids}")
print(f"Reaction Mode Symbol: {reaction_3.reaction_mode_symbol}")
print(
    f"reaction type: {reaction_3.reaction_type}"
)
print(
    f"Symbolic Unbalanced Reaction: {reaction_3.symbolic_unbalanced_reaction}")
print(f"Symbolic Reaction: {reaction_3.symbolic_reaction}")
print(f"Reactants: {reaction_3.reactants}")
print(f"Products: {reaction_3.products}")
print(f"Reaction Coefficients: {reaction_3.reaction_coefficients}")
print(f"Reaction Stoichiometry: {reaction_3.reaction_stoichiometry}")
print(
    f"Reaction Stoichiometry Matrix: {reaction_3.reaction_stoichiometry_matrix}")
print(
    f"Reaction Stoichiometry Source: {reaction_3.reaction_stoichiometry_source}")
print(f"State Counts: {reaction_3.state_count}")
print(f"Reaction Phase: {reaction_3.reaction_phase}")
print(f"Reaction State: {reaction_3.reaction_state}")
print(f"Carbon Count: {reaction_3.carbon_count}")
print(f"Reactants Names: {reaction_3.reactants_names}")
print(f"Products Names: {reaction_3.products_names}")
print(f"All Components: {reaction_3.all_components}")
print(f"Available components: {reaction_3.available_components}")
print(f"Component Checker: {reaction_3.component_checker}")
print(f"Mapped Components: {reaction_3.map_components}")
