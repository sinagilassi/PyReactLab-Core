# 📚 Documentation

Welcome to the PyReactLab-Core documentation! Here you'll find all the information you need to get started with PyReactLab-Core, including installation instructions, usage guides, and API references.

## 🧩 Package classes at a glance

The package is organized around a few core classes that capture chemical reactions, their structure, and the tools used to analyze them.

### `Reaction`

The [Reaction](reaction.md) model represents a single chemical equation and performs immediate validation and analysis. It parses reactants and products, derives stoichiometric information, checks element and charge balance, and exposes phase and component mapping.

Main API surface:

- `Reaction(name: str, reaction: str, balance_reaction: bool = False, components=None, reaction_mode_symbol=None, analysis={}, component_keys=[...])`
  - Parameters: `name` identifies the reaction, `reaction` is the equation string, `balance_reaction` enables automatic balancing, `components` optionally binds explicit component metadata, and `reaction_mode_symbol` sets the separator style.
  - Source: `pyreactlab_core.models.reaction.Reaction`

- `model_validate(...)` / model initialization validation
  - Notes: runs on instantiation and populates `self.analysis` before the object is used.

- Derived properties such as `reaction_stoichiometry`, `reaction_stoichiometry_matrix`, `is_element_balanced`, `is_charge_balanced`, and `is_balanced`
  - Notes: these are computed from the analyzed reaction payload and are the fields users inspect most often.

### `ReactionNetwork`

The [ReactionNetwork](reaction_network.md) model works with a set of reactions together. It builds the stoichiometric matrix, identifies independent versus dependent reactions, computes rank and dependency information, and evaluates network-level balance and grouping properties.

Main API surface:

- `ReactionNetwork(name: str, reactions: list[Reaction], balance_reaction: bool = False, components=None)`
  - Parameters: `name` is the network label, `reactions` is the set of reactions to analyze, `balance_reaction` optionally balances each reaction before the network is built, and `components` can fix the row ordering for the matrix.
  - Source: `pyreactlab_core.models.reaction_network.ReactionNetwork`

- `stoichiometric_matrix_by_components(components=None) -> list[list[float]]`
  - Notes: returns a matrix ordered by the supplied component list or by the network default order.

- `stoichiometric_matrix_dict_by_components(components=None) -> dict[str, list[float]]`
  - Notes: same as above, but keyed by component ID.

- `reaction_dependencies -> dict[str, dict[str, float]]`
  - Notes: expresses dependency relationships between reactions using the stoichiometric matrix.

### `ChemReact`

`ChemReact` is the core reaction-analysis engine behind the higher-level models. It performs parsing, normalization, phase/state inspection, balance-related checks, and structural interpretation of chemical equations.

Main methods and signatures:

- `analyze_reaction(reaction_pack: Dict[str, str], phase_rule: Optional[str] = None) -> Dict[str, Any]`
  - Parameters: `reaction_pack` must contain `name` and `reaction`; `phase_rule` optionally constrains the allowed state for the reaction.
  - Returns: a dictionary containing reactants, products, stoichiometry, phase, carbon counts, charge counts, and balance flags.
  - Source: `pyreactlab_core.core.chem_react.ChemReact.analyze_reaction`

- `parse_molecule(id: str, charge: Any) -> str`
  - Notes: formats a species as a molecular identifier, appending charge notation such as `{+}` or `{2-}` when needed.

- `count_carbon(molecule: str, coefficient: float) -> float`
  - Notes: counts carbon atoms in a species and scales the result by the reaction coefficient.

- `count_total_carbon(reactants: List[Reactant], products: List[Product]) -> Dict[str, float]`
  - Notes: returns total reactant carbon, total product carbon, and net carbon count.

- `phase_rule_analysis(phase_rule: Optional[str] = None) -> str`
  - Notes: normalizes the state rule (`gas`, `liquid`, `aqueous`, `solid`) into the compact reaction-state format the parser uses.

- `determine_reaction_phase(reaction_dict: Dict[str, str]) -> str`
  - Notes: resolves whether a reaction is single-phase or multiphase based on the species states.

- `get_reaction_type(reaction_mode_symbol: str) -> str`
  - Notes: classifies the equation as `reversible`, `irreversible`, or `equilibrium`.

- `parse_charge(charge: str) -> int`
  - Notes: converts charge notation to an integer value for the internal balance logic.

### `ChemReactBalance`

`ChemReactBalance` focuses on balancing reactions and related algebraic operations. It supports the conversion of unbalanced equations into a consistent stoichiometric form for subsequent analysis.

Main methods and signatures:

- `parse_ionic_charge(charge: str) -> int`
  - Notes: converts ionic charge strings into integer values used in element and charge bookkeeping.

- `parse_elemental_composition(formula: str) -> Dict[str, int]`
  - Notes: decomposes a molecule into elemental counts, including grouped ions and nested formula sections.

- `count_elements(molecule: str) -> Dict[str, int]`
  - Notes: counts the number of each element in a molecule formula.

- `elemental_balance(reactants: List[Reactant], products: List[Product]) -> Dict[str, Any]`
  - Notes: validates whether each element is conserved between the two sides of the reaction.

- `count_side_charge(side: List[Reactant] | List[Product]) -> Dict[str, Any]`
  - Notes: sums the charge contribution from each side of the reaction.

- `charge_balance(reactants: List[Reactant], products: List[Product]) -> Dict[str, Any]`
  - Notes: determines whether electrical charge is conserved in the reaction.

- `reaction_balance(reactants: List[Reactant], products: List[Product]) -> Dict[str, Any]`
  - Notes: consolidates elemental and charge checks into one final balance result.
  - Source: `pyreactlab_core.core.chem_react_balance.ChemReactBalance.reaction_balance`

### `ChemReactUtils`

`ChemReactUtils` contains reusable utilities for reaction parsing and data preparation, supporting cleaner handling of component identity, expression normalization, and analytical workflows.

Main methods and signatures:

- `parse_molecule(id: str, charge: Any) -> str`
  - Notes: builds a formatted molecule identifier including charge decoration when present.

- `count_carbon(molecule: str, coefficient: float) -> float`
  - Notes: counts carbon atoms in a formula and scales by the coefficient.

- `count_total_carbon(reactants, products) -> Dict[str, float]`
  - Notes: totals carbon on each side and reports net carbon transfer.

- `phase_rule_analysis(phase_rule: Optional[str] = None) -> str`
  - Notes: normalizes phase constraints for parsed species states.

- `state_name_set(state_set: set) -> List[str]`
  - Notes: converts compact reaction states to full names like `gas` and `aqueous`.

- `determine_reaction_phase(reaction_dict: Dict[str, str]) -> str`
  - Notes: infers the overall reaction phase from the set of species states.

- `count_reaction_states(reaction_dict: Dict[str, str]) -> Dict[str, int]`
  - Notes: counts occurrences of `g`, `l`, `aq`, and `s` in the relevant reaction species.

### `ReactionComponentMapper`

`ReactionComponentMapper` handles the mapping between reaction species and component definitions, helping connect parsed species to explicit component metadata when it is provided.

Main methods and signatures:

- `collect_components(reactants: Sequence[Any], products: Sequence[Any]) -> List[Component]`
  - Notes: gathers configured component objects that appear in the reaction and filters duplicates.

- `map_components(reactants: Sequence[Any], products: Sequence[Any]) -> Dict[str, Component]`
  - Notes: returns a lookup of `molecule_state` values to their configured `Component` objects.

- `reaction_stoichiometry_dict(reaction_stoichiometry: Dict[str, float], component_key: ComponentKey = "Name-Formula") -> Dict[str, float]`
  - Notes: rekeys stoichiometry values using a requested metadata key such as `Name-Formula`, `Formula-State`, or similar.

- `build_stoichiometry_source(reaction_stoichiometry: Dict[str, float]) -> Dict[str, Any]`
  - Notes: builds the full component-key mapping used by the model for downstream analysis.
  - Source: `pyreactlab_core.core.reaction_component_mapper.ReactionComponentMapper.build_stoichiometry_source`

## 🔗 Related pages

- [Reaction Model](reaction.md)
- [Reaction Network Model](reaction_network.md)
- [Models Overview](models.md)


## Reaction Analysis

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

## Reaction Network Analysis

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

## ⚗️ Reaction Analysis made easy

::: pyreactlab_core.core.chem_react
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

## ⚖️ Reaction Stoichiometry made easy

::: pyreactlab_core.docs.chem_balance
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