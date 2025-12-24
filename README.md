# PyReactLab-Core

[![PyPI Downloads](https://static.pepy.tech/badge/pyreactlab-core/month)](https://pepy.tech/projects/pyreactlab-core)
![PyPI](https://img.shields.io/pypi/v/pyreactlab-core)
![Python Version](https://img.shields.io/pypi/pyversions/pyreactlab-core.svg)
![License](https://img.shields.io/pypi/l/pyreactlab-core)

**PyReactLab-Core** is the core foundation of the PyReactLab ecosystem, offering shared data structures and algorithms for chemical reaction representation, stoichiometry, and reaction analysis.

## Features

- **Reaction Representation**: Define and manipulate chemical reactions with ease.
- **Stoichiometry Calculations**: Perform stoichiometric calculations for reactions.
- **Reaction Analysis**: Analyze reaction properties and behaviors.
- **Extensible Design**: Built to be extended by other PyReactLab modules.

## Installation

You can install PyReactLab-Core via pip:

```bash
pip install pyreactlab-core
```

## 🚀 Usage

### Introduce a reaction

A typical reaction can be introduced as follows:

```python
from pyreactlab_core import Reaction

reaction = Reaction(
    name="Combustion of Methane",
    reaction="CO2(g) + 3H2(g) => CH3OH(g) + H2O(g)"
)

# print analysis
print(
    f"[bold underline]Reaction Analysis for: {reaction_1.name}[/bold underline]")
print(f"Reaction: {reaction_1.reaction}")
print(f"Reactants: {reaction_1.reactants}")
print(f"Products: {reaction_1.products}")
print(f"Reaction Coefficients: {reaction_1.reaction_coefficients}")
print(f"Reaction Stoichiometry: {reaction_1.reaction_stoichiometry}")
print(f"State Counts: {reaction_1.state_count}")
print(f"Reaction Phase: {reaction_1.reaction_phase}")
print(f"Reaction State: {reaction_1.reaction_state}")
print(f"Carbon Count: {reaction_1.carbon_count}")
print(f"Reactants Names: {reaction_1.reactants_names}")
print(f"Products Names: {reaction_1.products_names}")

# results:
# Reaction: CO2(g) + 3H2(g) => CH3OH(g) + H2O(g)
# Component IDs: {'CO2-g': 1, 'H2-g': 2, 'CH3OH-g': 3, 'H2O-g': 4}
# Reaction Mode Symbol: =>
# Symbolic Unbalanced Reaction: CO2 + H2 => CH3OH + H2O
# Symbolic Reaction: CO2 + 3.0H2 => CH3OH + H2O
# Reactants: [{'coefficient': 1.0, 'molecule': 'CO2', 'state': 'g', 'molecule_state': 'CO2-g'}, {'coefficient': 3.0, 'molecule': 'H2', 'state': 'g', 'molecule_state': 'H2-g'}]
# Products: [{'coefficient': 1.0, 'molecule': 'CH3OH', 'state': 'g', 'molecule_state': 'CH3OH-g'}, {'coefficient': 1.0, 'molecule': 'H2O', 'state': 'g', 'molecule_state': 'H2O-g'}]
# Reaction Coefficients: 2.0
# Reaction Stoichiometry: {'CO2-g': -1.0, 'H2-g': -3.0, 'CH3OH-g': 1.0, 'H2O-g': 1.0}
# State Counts: {'g': 4, 'l': 0, 'aq': 0, 's': 0}
# Reaction Phase: gas
# Reaction State: {'CO2-g': 'g', 'H2-g': 'g', 'CH3OH-g': 'g', 'H2O-g': 'g'}
# Carbon Count: {'CO2-g': 1.0, 'H2-g': 0.0, 'CH3OH-g': 1.0, 'H2O-g': 0.0}
# Reactants Names: ['CO2-g', 'H2-g']
# Products Names: ['CH3OH-g', 'H2O-g']
```

### Stoichiometric Balance

You can check if a reaction is balanced:

```python
from pyreactlab_core import Reaction
from pyreactlab_core.core import balance

# define a reaction
reaction = Reaction(
    name="Combustion of Methane",
    reaction="CO2(g) + 3H2(g) => CH3OH(g) + H2O(g)"
)

# balance the reaction automatically
balanced_reaction = balance(reaction)
print(f"Balanced Reaction: {balanced_reaction.reaction}")
```

### Stoichiometric Matrix

You can create a stoichiometric matrix for a list of reactions:

```python
from pyreactlab_core import Reaction
from pyreactlab_core import rxn, rxn_stoichiometry, rxns_stoichiometry

# NOTE: define reaction string
reaction_1 = "CO2(g) + 3H2(g) => CH3OH(g) + H2O(g)"
name_1 = "CO2 Hydrogenation to Methanol"

# second reaction
reaction_2 = "C2H4(g) + H2(g) => C2H6(g)"
name_2 = "Ethylene Hydrogenation to Ethane"

# NOTE: create reaction instance
rxn_1: Reaction = rxn(
    reaction_str=reaction_1,
    name=name_1
)

rxn_2: Reaction = rxn(
    reaction_str=reaction_2,
    name=name_2
)

# NOTE: Get stoichiometry matrices for multiple reactions
reactions_list = [rxn_1, rxn_2]
stoichiometry_matrices = rxns_stoichiometry(
    reactions=reactions_list,
)
# log
print(stoichiometry_matrices)

# results:
# {
#     'components': ['CO2-g', 'H2O-g', 'C2H4-g', 'H2-g', 'CH3OH-g', 'C2H6-g'],
#     'component_ids': {'CO2-g': 0, 'H2O-g': 1, 'C2H4-g': 2, 'H2-g': 3, 'CH3OH-g': 4, 'C2H6-g': 5},
#     'stoichiometry_matrices_list': [[-1.0, 1.0, 0.0, -3.0, 1.0, 0.0], [0.0, 0.0, -1.0, -1.0, 0.0, 1.0]],
#     'stoichiometry_matrices_dict': [
#         {'CO2-g': -1.0, 'H2O-g': 1.0, 'C2H4-g': 0.0, 'H2-g': -3.0, 'CH3OH-g': 1.0, 'C2H6-g': 0.0},
#         {'CO2-g': 0.0, 'H2O-g': 0.0, 'C2H4-g': -1.0, 'H2-g': -1.0, 'CH3OH-g': 0.0, 'C2H6-g': 1.0}
#     ]
# }
```

## 🤝 Contributing

Contributions are highly welcome — bug fixes, new calculation routines, mixture models, extended unit tests, documentation, etc.

## 📝 License

This project is distributed under the Apache License, Version 2.0, which grants you broad freedom to use, modify, and integrate the software into your own applications or projects, provided that you comply with the conditions outlined in the license. Although Apache 2.0 does not require users to retain explicit author credit beyond standard copyright and license notices, I kindly request that if you incorporate this work into your own software, you acknowledge Sina Gilassi as the original author. Referencing the original repository or documentation is appreciated, as it helps recognize the effort invested in developing and maintaining this project.

## ❓ FAQ

For any question, contact me on [LinkedIn](https://www.linkedin.com/in/sina-gilassi/)

## 👨‍💻 Authors

- [@sinagilassi](https://www.github.com/sinagilassi)
