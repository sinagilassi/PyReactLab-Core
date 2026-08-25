# import libs
from pythermodb_settings.models import Component
from pyreactlab_core.models import Reaction, ReactionNetwork
import sys
from pathlib import Path
from rich import print

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))


# SECTION: define components
component_h2 = Component(
    name="Hydrogen",
    formula="H2",
    state="g",
)

component_h = Component(
    name="Atomic Hydrogen",
    formula="H",
    state="g",
)

component_h2o = Component(
    name="Water",
    formula="H2O",
    state="g",
)

component_o = Component(
    name="Atomic Oxygen",
    formula="O",
    state="g",
)

components = [
    component_h2,
    component_h,
    component_h2o,
    component_o,
]


# SECTION: define reactions
# NOTE: H-g is produced by R1 and consumed by R2, so it is an intermediate.
reaction_1 = Reaction(
    name="R1",
    reaction="H2(g) => 2H(g)",
    components=[
        component_h2,
        component_h,
    ],
)

reaction_2 = Reaction(
    name="R2",
    reaction="2H(g) + O(g) => H2O(g)",
    components=[
        component_h,
        component_o,
        component_h2o,
    ],
)


# SECTION: build reaction network
network = ReactionNetwork(
    name="intermediate-species-network",
    reactions=[
        reaction_1,
        reaction_2,
    ],
)


# SECTION: display network identity
print("Intermediate Species Network:")
print(f"  reaction_ids: {network.reaction_ids}")
print(f"  component_ids: {network.component_ids}")
print(
    f"  configured_components: {[f'{item.formula}-{item.state}' for item in components]}")


# SECTION: display stoichiometry
print("Stoichiometric Matrix:")
for component_id, row in zip(network.component_ids, network.stoichiometric_matrix):
    print(f"  {component_id}: {row}")


# SECTION: display species role classification
print("Species Roles:")
print(f"  reactants: {network.reactants}")
print(f"  products: {network.products}")
print(f"  source_species: {network.source_species}")
print(f"  intermediate_species: {network.intermediate_species}")
print(f"  sink_species: {network.sink_species}")


# SECTION: display structural analysis
print("Structural Analysis:")
print(f"  stoichiometric_rank: {network.stoichiometric_rank}")
print(f"  independent_reactions: {network.independent_reactions}")
print(f"  dependent_reactions: {network.dependent_reactions}")
print(f"  reaction_dependencies: {network.reaction_dependencies}")


# SECTION: display summary
print("Summary:")
print(network.summary)
