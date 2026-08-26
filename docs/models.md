# 🧬 Models

This section provides a high-level overview of the core data models used in PyReactLab-Core. The project centers on two principal model types:

- [Reaction](reaction.md): a single chemical reaction with parsed stoichiometry, phase metadata, and balance analysis
- [ReactionNetwork](reaction_network.md): a collection of reactions analyzed together as a stoichiometric network

## 🔗 Related pages

- [Reaction Model](reaction.md)
- [Reaction Network Model](reaction_network.md)
- [Documentation Home](docs.md)

## Overview

`Reaction` represents one equation and exposes its structural analysis immediately after initialization, including reactant/product parsing, stoichiometric coefficients, element and charge balance checks, and component mapping.

`ReactionNetwork` builds on this by analyzing multiple `Reaction` objects together. It derives the stoichiometric matrix, detects independent and dependent reactions, groups species by role and phase, and evaluates element and charge balance across the full network.

## Model relationships

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

network = ReactionNetwork(
    name="methanol-synthesis",
    reactions=[r1, r2],
)

print(r1.is_balanced)
print(network.stoichiometric_rank)
print(network.independent_reactions)
```

The detailed API documentation for each model lives on its dedicated page to avoid repeating the same reference material in the overview page.
