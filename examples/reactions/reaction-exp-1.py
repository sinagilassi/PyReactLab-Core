# import libs
from pyreactlab_core.models.reaction import Reaction
from pythermodb_settings.models import Component
from rich import print

# NOTE: define components
component_co2 = Component(
    name="Carbon Dioxide",
    formula="CO2",
    state="g"
)

component_h2 = Component(
    name="Hydrogen",
    formula="H2",
    state="g"
)

component_ch3oh = Component(
    name="Methanol",
    formula="CH3OH",
    state="g"
)

component_h2o = Component(
    name="Water",
    formula="H2O",
    state="g"
)

component_co = Component(
    name="Carbon Monoxide",
    formula="CO",
    state="g"
)

components = [
    component_h2,
    component_co,
    component_ch3oh,
    component_h2o,
    component_co2,
]

# components = []
# components = None

# NOTE: define components
reaction_1 = Reaction(
    name="Combustion of Methane",
    reaction="CO2(g) + 3H2(g) ⇌ CH3OH(g) + H2O(g)",
    components=components
)

# NOTE: print analysis
print(
    f"[bold underline]Reaction Analysis for: {reaction_1.name}[/bold underline]")
print("Reaction: ")
print(reaction_1.reaction)
print("Component IDs: ")
print(reaction_1.component_ids)
print("Reaction Mode Symbol: ")
print(reaction_1.reaction_mode_symbol)
print("Reaction Type: ")
print(reaction_1.reaction_type)
print("Symbolic Unbalanced Reaction: ")
print(reaction_1.symbolic_unbalanced_reaction)
print("Symbolic Reaction: ")
print(reaction_1.symbolic_reaction)
print("Reactants: ")
print(reaction_1.reactants)
print("Products: ")
print(reaction_1.products)
print("Reaction Coefficients: ")
print(reaction_1.reaction_coefficients)
print("Reaction Stoichiometry: ")
print(reaction_1.reaction_stoichiometry)
print("Reaction Stoichiometry Matrix: ")
print(reaction_1.reaction_stoichiometry_matrix)
print("Reaction Stoichiometry Source: ")
print(reaction_1.reaction_stoichiometry_source)
print("State Counts: ")
print(reaction_1.state_count)
print("Reaction Phase: ")
print(reaction_1.reaction_phase)
print("Reaction State: ")
print(reaction_1.reaction_state)
print("Carbon Count: ")
print(reaction_1.carbon_count)
print("Reactants Names: ")
print(reaction_1.reactants_names)
print("Products Names: ")
print(reaction_1.products_names)
print("All Components: ")
print(reaction_1.all_components)
print("Available components: ")
print(reaction_1.available_components)
print("Component Checker: ")
print(reaction_1.component_checker)
print("Mapped Components: ")
print(reaction_1.map_components)
print(f"Charge Count: ")
print(reaction_1.charge_count)
print("total reactant charge: ")
print(reaction_1.total_reactant_charge)
print("total product charge: ")
print(reaction_1.total_product_charge)
print("net charge: ")
print(reaction_1.net_charge)
