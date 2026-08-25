# import libs
from pythermodb_settings.models import Component
from pyreactlab_core.models import Reaction, ReactionNetwork
import sys
from pathlib import Path
from rich import print

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


# SECTION: build reaction network
network = ReactionNetwork(
    name="methanol-synthesis-network",
    reactions=[
        reaction_1,
        reaction_2,
        reaction_3,
    ],
)


# SECTION: display network identity and stoichiometry
print("Methanol Synthesis Network:")
print(f"  reaction_ids: {network.reaction_ids}")
print(f"  component_ids: {network.component_ids}")
print(
    f"  configured_components: {[f'{item.formula}-{item.state}' for item in components]}")
print("  stoichiometric_matrix:")
for component_id, row in zip(network.component_ids, network.stoichiometric_matrix):
    print(f"    {component_id}: {row}")


# SECTION: display dependency analysis
# NOTE: The result should show R3 = R1 - R2.
print("Dependency Analysis:")
print(f"  stoichiometric_rank: {network.stoichiometric_rank}")
print(f"  independent_reactions: {network.independent_reactions}")
print(f"  dependent_reactions: {network.dependent_reactions}")
print(f"  reaction_dependencies: {network.reaction_dependencies}")


# SECTION: display structural maps
print("Structural Maps:")
print(f"  reaction_component_map: {network.reaction_component_map}")
print(f"  component_reaction_map: {network.component_reaction_map}")
print(f"  participation_matrix: {network.participation_matrix}")


# SECTION: display phase and balance data
print("Phase and Balance:")
print(f"  phases: {network.phases}")
print(f"  components_by_phase: {network.components_by_phase}")
print(f"  is_element_balanced: {network.is_element_balanced}")
print(f"  is_charge_balanced: {network.is_charge_balanced}")
print(f"  is_balanced: {network.is_balanced}")


# SECTION: display summary
print("Summary:")
print(network.summary)
