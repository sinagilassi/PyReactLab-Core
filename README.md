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

## 🤝 Contributing

Contributions are highly welcome — bug fixes, new calculation routines, mixture models, extended unit tests, documentation, etc.

## 📝 License

This project is distributed under the Apache License, Version 2.0, which grants you broad freedom to use, modify, and integrate the software into your own applications or projects, provided that you comply with the conditions outlined in the license. Although Apache 2.0 does not require users to retain explicit author credit beyond standard copyright and license notices, I kindly request that if you incorporate this work into your own software, you acknowledge Sina Gilassi as the original author. Referencing the original repository or documentation is appreciated, as it helps recognize the effort invested in developing and maintaining this project.

## ❓ FAQ

For any question, contact me on [LinkedIn](https://www.linkedin.com/in/sina-gilassi/)

## 👨‍💻 Authors

- [@sinagilassi](https://www.github.com/sinagilassi)