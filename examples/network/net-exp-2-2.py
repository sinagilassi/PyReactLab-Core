# import libs
from rich import print
from pyreactlab_core.models import Reaction, ReactionNetwork
from pythermodb_settings.models import Component
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))


# SECTION: define components
component_co2 = Component(
    name="Carbon Dioxide",
    formula="CO2",
    state="g",
)

component_h2 = Component(
    name="Hydrogen",
    formula="H2",
    state="g",
)

component_co = Component(
    name="Carbon Monoxide",
    formula="CO",
    state="g",
)

component_ch3oh = Component(
    name="Methanol",
    formula="CH3OH",
    state="g",
)

component_h2o = Component(
    name="Water",
    formula="H2O",
    state="g",
)

components = [
    component_co2,
    component_h2,
    component_co,
    component_ch3oh,
    component_h2o,
]

reordered_components = [
    component_h2o,
    component_ch3oh,
    component_co,
    component_h2,
    component_co2,
]


# SECTION: define methanol synthesis network reactions
# NOTE: With this ordering, R1 and R2 form the selected basis.
reaction_1 = Reaction(
    name="R1",
    reaction="CO2(g) + 3H2(g) <=> CH3OH(g) + H2O(g)",
    components=[
        component_co2,
        component_h2,
        component_ch3oh,
        component_h2o,
    ],
)

reaction_2 = Reaction(
    name="R2",
    reaction="CO2(g) + H2(g) <=> CO(g) + H2O(g)",
    components=[
        component_co2,
        component_h2,
        component_co,
        component_h2o,
    ],
)

reaction_3 = Reaction(
    name="R3",
    reaction="CO(g) + 2H2(g) <=> CH3OH(g)",
    components=[
        component_co,
        component_h2,
        component_ch3oh,
    ],
)


# SECTION: build reaction networks
network_default = ReactionNetwork(
    name="methanol-synthesis-network-default",
    reactions=[
        reaction_1,
        reaction_2,
        reaction_3,
    ],
)

network_reordered = ReactionNetwork(
    name="methanol-synthesis-network",
    reactions=[
        reaction_1,
        reaction_2,
        reaction_3,
    ],
    components=reordered_components,
)


# SECTION: compare default and configured component order
print("Methanol Synthesis Network:")
print(f"  reaction_ids: {network_default.reaction_ids}")
print(
    f"  configured_components: {[f'{item.formula}-{item.state}' for item in components]}")
print(
    "  reordered_components: "
    f"{[f'{item.formula}-{item.state}' for item in reordered_components]}"
)

print("Default Network (components=None):")
print(f"  component_ids: {network_default.component_ids}")
print("  stoichiometric_matrix:")
for component_id, row in zip(
    network_default.component_ids,
    network_default.stoichiometric_matrix,
):
    print(f"    {component_id}: {row}")
print("  stoichiometric_matrix_dict:")
print(network_default.stoichiometric_matrix_dict)
print(" stoichiometric_matrix:")
print(network_default.stoichiometric_matrix)

print("Configured Network (components=reordered_components):")
print(f"  component_ids: {network_reordered.component_ids}")
print("  stoichiometric_matrix:")
for component_id, row in zip(
    network_reordered.component_ids,
    network_reordered.stoichiometric_matrix,
):
    print(f"    {component_id}: {row}")
print("  stoichiometric_matrix_dict:")
print(network_reordered.stoichiometric_matrix_dict)
print("  stoichiometric matrix:")
print(network_reordered.stoichiometric_matrix)
