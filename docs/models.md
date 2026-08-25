# 🧬 Models

This section provides an overview of the core data models used in PyReactLab-Core, including classes for representing chemical species, reactions, and reaction networks. Each model is designed to facilitate efficient storage, manipulation, and analysis of chemical data.

## ⚛️ Chemical Reaction Model

::: pyreactlab_core.models.reaction
    handler: python
    options:
        show_source: true
        show_root_heading: true
        show_root_full_path: true
        show_signature_annotations: true
        docstring_style: google
        merge_init_into_class: true
        show_if_no_docstring: true
        show_bases: true
        show_inheritance_diagram: false
        show_submodules: true

## Reaction Network Model

`ReactionNetwork` represents and analyzes relationships among multiple
`Reaction` objects. It stays within structural reaction analysis: species
collection, stoichiometry, rank, reaction dependency, phases, charge, and
elemental consistency. Operating conditions, reaction extents, equilibrium,
kinetics, reactor models, optimization, and simulation belong downstream in
PyReactLab's `ReactionSystem`.

The stoichiometric matrix convention is explicit:

- rows = components/species
- columns = reactions
- reactant coefficients are negative
- product coefficients are positive
- non-participating species are zero

Component identity uses the parsed molecule-state identifier exposed by
`Reaction`, such as `CO2-g`, `H2O-l`, `H{+}-aq`, `OH{-}-aq`, and
`Fe{3+}-aq`. Formula, charge, and phase are therefore preserved in network
matrices and maps.

```python
from pyreactlab_core.models import Reaction, ReactionNetwork

r1 = Reaction(
    name="R1",
    reaction="CO2(g) + 3H2(g) <=> CH3OH(g) + H2O(g)",
)
r2 = Reaction(
    name="R2",
    reaction="CO2(g) + H2(g) <=> CO(g) + H2O(g)",
)
r3 = Reaction(
    name="R3",
    reaction="CO(g) + 2H2(g) <=> CH3OH(g)",
)

network = ReactionNetwork(
    name="methanol-synthesis",
    reactions=[r1, r2, r3],
)

print(network.component_ids)
print(network.stoichiometric_matrix)
print(network.stoichiometric_rank)
print(network.independent_reactions)
print(network.dependent_reactions)
print(network.reaction_dependencies)
print(network.summary)
```

::: pyreactlab_core.models.reaction_network
    handler: python
    options:
        show_source: true
        show_root_heading: true
        show_root_full_path: true
        show_signature_annotations: true
        docstring_style: google
        merge_init_into_class: true
        show_if_no_docstring: true
        show_bases: true
        show_inheritance_diagram: false
        show_submodules: true
