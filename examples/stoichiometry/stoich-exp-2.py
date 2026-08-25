# import libs
from pyreactlab_core.models.reaction import Reaction
from pythermodb_settings.models import Component
from rich import print
from pyreactlab_core import build_rxns_stoichiometry

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

component_c2h4 = Component(
    name="Ethylene",
    formula="C2H4",
    state="g"
)

component_c2h6 = Component(
    name="Ethane",
    formula="C2H6",
    state="g"
)

component_co = Component(
    name="Carbon Monoxide",
    formula="CO",
    state="g"
)

components = [
    component_co2,
    component_h2,
    component_ch3oh,
    component_co,
    component_h2o,
    component_c2h4,
    component_c2h6
]

# NOTE: define reaction string
reaction_1 = "CO2(g) + 3H2(g) => CH3OH(g) + H2O(g)"
name_1 = "CO2 Hydrogenation to Methanol"
components_1 = [component_co2, component_h2, component_ch3oh, component_h2o]

# second reaction
reaction_2 = "C2H4(g) + H2(g) => C2H6(g)"
name_2 = "Ethylene Hydrogenation to Ethane"
components_2 = [component_c2h4, component_h2, component_c2h6]

# NOTE: define reactions
reaction_1 = Reaction(
    name=name_1,
    reaction=reaction_1,
    components=components_1
)

reaction_2 = Reaction(
    name=name_2,
    reaction=reaction_2,
    components=components_2
)

# SECTION: print analysis
# Reaction 1
print(
    f"[bold underline]Reaction Analysis for: {reaction_1.name}[/bold underline]")
print(f"Reaction: {reaction_1.reaction}")
print(f"Component IDs: {reaction_1.component_ids}")
print(f"Reaction Coefficients: {reaction_1.reaction_coefficients}")
print(f"Reaction Stoichiometry: {reaction_1.reaction_stoichiometry}")
print(
    f"Reaction Stoichiometry Matrix: {reaction_1.reaction_stoichiometry_matrix}")
print(
    f"Reaction Stoichiometry Source: {reaction_1.reaction_stoichiometry_source}")

# Reaction 2
print(
    f"\n[bold underline]Reaction Analysis for: {reaction_2.name}[/bold underline]")
print(f"Reaction: {reaction_2.reaction}")
print(f"Component IDs: {reaction_2.component_ids}")
print(f"Reaction Coefficients: {reaction_2.reaction_coefficients}")
print(f"Reaction Stoichiometry: {reaction_2.reaction_stoichiometry}")
print(
    f"Reaction Stoichiometry Matrix: {reaction_2.reaction_stoichiometry_matrix}")
print(
    f"Reaction Stoichiometry Source: {reaction_2.reaction_stoichiometry_source}")

# Components
print("\n[bold underline]Components List[/bold underline]")
print(components)


# SECTION: Get reaction stoichiometry matrix
stoichiometry_result = build_rxns_stoichiometry(
    reactions=[reaction_1, reaction_2],
    components=components,
    component_key="Name-Formula"
)
print(stoichiometry_result)
