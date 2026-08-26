# 🧩 Reaction Network Model

This page documents the `ReactionNetwork` model, which analyzes how multiple chemical reactions relate to one another structurally.

`ReactionNetwork` focuses on reaction topology and stoichiometric relationships: species participation, independent versus dependent reactions, matrix rank, phase grouping, and global element/charge consistency. It is the structural core used before downstream reactor, equilibrium, or kinetic modeling.

The stoichiometric matrix convention is explicit:

- rows = components/species
- columns = reactions
- reactant coefficients are negative
- product coefficients are positive
- non-participating species are zero

Component identity uses the parsed molecule-state identifier exposed by `Reaction`, such as `CO2-g`, `H2O-l`, `H{+}-aq`, `OH{-}-aq`, and `Fe{3+}-aq`. Formula, charge, and phase are therefore preserved in network matrices and maps.

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

## Minimal Example

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

## Core Concepts

`ReactionNetwork` provides:

- species classification (source, sink, intermediate)
- stoichiometric matrix and rank
- independent versus dependent reaction detection
- reaction dependency coefficients
- phase and reaction mode grouping
- element and charge balance checks

## Interpreting Dependency Results

If the dependency output looks like:

```python
{"R3": {"R1": 1.0, "R2": -1.0}}
```

it means that `R3` is dependent and can be expressed mathematically as:

$$
R3 = 1.0 \times R1 - 1.0 \times R2
$$

This is derived from the column-space linear dependency of the network stoichiometric matrix.

## Network Properties You Will Use Most

The most commonly consulted properties are:

- `component_ids`
- `stoichiometric_matrix`
- `stoichiometric_rank`
- `independent_reactions`
- `dependent_reactions`
- `reaction_dependencies`
- `summary`

These give a compact but powerful view of the structural relationships among reactions in a network.
