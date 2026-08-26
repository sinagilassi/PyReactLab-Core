# ⚗️ Reaction Model

This page documents the core `Reaction` model used to represent a single chemical transformation, parse its stoichiometry, and evaluate its structural balance properties.

The model is designed to work directly with reaction strings such as `CO2(g) + 3H2(g) <=> CH3OH(g) + H2O(g)` and then expose parsed reactant/product data, elemental balance checks, charge analysis, and component mapping in a consistent format.

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

## Overview

`Reaction` is a Pydantic model that captures a single chemical equation and analyzes it immediately at initialization time.

It exposes:

- normalized reaction text
- reactant and product parsing
- stoichiometric coefficients and matrix data
- phase/state information
- element and charge balance checks
- component mapping and validation

## Typical Usage

```python
from pyreactlab_core.models import Reaction

reaction = Reaction(
    name="R1",
    reaction="CO2(g) + 3H2(g) <=> CH3OH(g) + H2O(g)",
)

print(reaction.symbolic_reaction)
print(reaction.reaction_stoichiometry)
print(reaction.is_element_balanced)
print(reaction.is_charge_balanced)
print(reaction.is_balanced)
```

## What to Read First

Most users start with these fields:

- `symbolic_reaction`
- `reaction_stoichiometry`
- `reaction_stoichiometry_matrix`
- `is_element_balanced`
- `is_charge_balanced`
- `is_balanced`

If you need the complete structured analysis payload, use:

```python
payload = reaction.model_dump()
```
