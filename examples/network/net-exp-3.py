# import libs
from pythermodb_settings.models import Component
from pyreactlab_core.models import Reaction, ReactionNetwork
import sys
from pathlib import Path
from rich import print

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))


# SECTION: define ionic components
component_h_plus = Component(
    name="Hydrogen Ion",
    formula="H{+}",
    state="aq",
)

component_oh_minus = Component(
    name="Hydroxide Ion",
    formula="OH{-}",
    state="aq",
)

component_h2o_l = Component(
    name="Water",
    formula="H2O",
    state="l",
)

component_fe3 = Component(
    name="Iron(III) Ion",
    formula="Fe{3+}",
    state="aq",
)

component_fe2 = Component(
    name="Iron(II) Ion",
    formula="Fe{2+}",
    state="aq",
)

component_electron = Component(
    name="Electron",
    formula="e{-}",
    state="aq",
)

components = [
    component_h_plus,
    component_oh_minus,
    component_h2o_l,
    component_fe3,
    component_fe2,
    component_electron,
]


# SECTION: define ionic reactions
neutralization = Reaction(
    name="R1",
    reaction="H{+}(aq) + OH{-}(aq) => H2O(l)",
    components=[
        component_h_plus,
        component_oh_minus,
        component_h2o_l,
    ],
)

iron_reduction = Reaction(
    name="R2",
    reaction="Fe{3+}(aq) + e{-}(aq) => Fe{2+}(aq)",
    components=[
        component_fe3,
        component_electron,
        component_fe2,
    ],
)


# SECTION: build ionic reaction network
network = ReactionNetwork(
    name="ionic-reaction-network",
    reactions=[
        neutralization,
        iron_reduction,
    ],
)


# SECTION: display ionic identity and stoichiometry
print("Ionic Reaction Network:")
print(f"  reaction_ids: {network.reaction_ids}")
print(f"  component_ids: {network.component_ids}")
print(
    f"  configured_components: {[f'{item.formula}-{item.state}' for item in components]}")
print("  stoichiometric_matrix:")
for component_id, row in zip(network.component_ids, network.stoichiometric_matrix):
    print(f"    {component_id}: {row}")


# SECTION: display charge-aware network data
print("Charge and Balance:")
print(f"  charge_vector: {network.charge_vector}")
print(f"  charge_balance_vector: {network.charge_balance_vector}")
print(f"  is_charge_balanced: {network.is_charge_balanced}")
print(f"  is_element_balanced: {network.is_element_balanced}")
print(f"  is_balanced: {network.is_balanced}")


# SECTION: display phase and dependency data
print("Phase and Dependency:")
print(f"  phases: {network.phases}")
print(f"  components_by_phase: {network.components_by_phase}")
print(f"  independent_reactions: {network.independent_reactions}")
print(f"  dependent_reactions: {network.dependent_reactions}")
print(f"  reaction_dependencies: {network.reaction_dependencies}")


# SECTION: display component mappings
print("Mapped Components:")
for reaction in network.reactions:
    print(f"  {reaction.name}: {list(reaction.map_components.keys())}")


# SECTION: display summary
print("Summary:")
print(network.summary)
