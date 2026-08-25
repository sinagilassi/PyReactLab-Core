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

component_o2 = Component(
    name="Oxygen",
    formula="O2",
    state="g",
)

component_h = Component(
    name="Atomic Hydrogen",
    formula="H",
    state="g",
)

component_o = Component(
    name="Atomic Oxygen",
    formula="O",
    state="g",
)

components = [
    component_h2,
    component_o2,
    component_h,
    component_o,
]


# SECTION: define reactions
# NOTE: R3 is stoichiometrically dependent because R3 = R1 + R2.
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
    reaction="O2(g) => 2O(g)",
    components=[
        component_o2,
        component_o,
    ],
)

reaction_3 = Reaction(
    name="R3",
    reaction="H2(g) + O2(g) => 2H(g) + 2O(g)",
    components=components,
)


# SECTION: build reaction network
network = ReactionNetwork(
    name="additive-dependency-network",
    reactions=[
        reaction_1,
        reaction_2,
        reaction_3,
    ],
)


# SECTION: display structural network data
print("Reaction Network:")
print(f"  name: {network.name}")
print(f"  reaction_ids: {network.reaction_ids}")
print(f"  component_ids: {network.component_ids}")
print(
    f"  configured_components: {[f'{item.formula}-{item.state}' for item in components]}")
print(f"  mapped_components_R3: {list(reaction_3.map_components.keys())}")
print(f"  stoichiometric_matrix_shape: {network.stoichiometric_matrix_shape}")
print("  stoichiometric_matrix:")
for row in network.stoichiometric_matrix:
    print(f"    {row}")


# SECTION: display dependency analysis
print("Dependency Analysis:")
print(f"  stoichiometric_rank: {network.stoichiometric_rank}")
print(
    f"  independent_reaction_indices: {network.independent_reaction_indices}")
print(f"  dependent_reaction_indices: {network.dependent_reaction_indices}")
print(f"  independent_reactions: {network.independent_reactions}")
print(f"  dependent_reactions: {network.dependent_reactions}")
print(
    f"  independent_stoichiometric_matrix: {network.independent_stoichiometric_matrix}")
print(f"  reaction_dependencies: {network.reaction_dependencies}")


# SECTION: display network species roles
print("Species Roles:")
print(f"  source_species: {network.source_species}")
print(f"  intermediate_species: {network.intermediate_species}")
print(f"  sink_species: {network.sink_species}")


# SECTION: display summary
print("Summary:")
print(network.summary)
