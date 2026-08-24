# import libs
from pyreactlab_core.models.reaction import Reaction
from pythermodb_settings.models import Component
from rich import print

# NOTE: C-starting element, but NOT carbon (Fe/Ce both change oxidation state,
# so the plain-formula identity used for charged species collides and no
# component list can uniquely map them)
reaction_1 = Reaction(
    name="Cerium-Iron Redox Reaction",
    reaction="Fe{2+}(aq) + Ce{4+}(aq) => Fe{3+}(aq) + Ce{3+}(aq)",
    components=None,
)

# NOTE: print analysis
print(
    f"[bold underline]Reaction Analysis for: {reaction_1.name}[/bold underline]")
print(f"Reaction: {reaction_1.reaction}")
print(f"Component IDs: {reaction_1.component_ids}")
print(f"Reaction Mode Symbol: {reaction_1.reaction_mode_symbol}")
print(
    f"reaction type: {reaction_1.reaction_type}"
)
print(
    f"Symbolic Unbalanced Reaction: {reaction_1.symbolic_unbalanced_reaction}")
print(f"Symbolic Reaction: {reaction_1.symbolic_reaction}")
print(f"Reactants: {reaction_1.reactants}")
print(f"Products: {reaction_1.products}")
print(f"Reaction Coefficients: {reaction_1.reaction_coefficients}")
print(f"Reaction Stoichiometry: {reaction_1.reaction_stoichiometry}")
print(
    f"Reaction Stoichiometry Matrix: {reaction_1.reaction_stoichiometry_matrix}")
print(
    f"Reaction Stoichiometry Source: {reaction_1.reaction_stoichiometry_source}")
print(f"State Counts: {reaction_1.state_count}")
print(f"Reaction Phase: {reaction_1.reaction_phase}")
print(f"Reaction State: {reaction_1.reaction_state}")
print(f"Carbon Count: {reaction_1.carbon_count}")
print(f"Reactants Names: {reaction_1.reactants_names}")
print(f"Products Names: {reaction_1.products_names}")
print(f"All Components: {reaction_1.all_components}")
print(f"Available components: {reaction_1.available_components}")
print(f"Component Checker: {reaction_1.component_checker}")
print(f"Mapped Components: {reaction_1.map_components}")

# NOTE: redox reaction (distinct states keep Zn/Cu identities unique)
component_zn_s = Component(
    name="Zinc",
    formula="Zn",
    state="s"
)
component_cu_ion = Component(
    name="Copper(II) Ion",
    formula="Cu",
    state="aq"
)
component_zn_ion = Component(
    name="Zinc Ion",
    formula="Zn",
    state="aq"
)
component_cu_s = Component(
    name="Copper",
    formula="Cu",
    state="s"
)

components_2 = [
    component_zn_s,
    component_cu_ion,
    component_zn_ion,
    component_cu_s,
]

reaction_2 = Reaction(
    name="Zinc-Copper Redox Reaction",
    reaction="Zn(s) + Cu{2+}(aq) => Zn{2+}(aq) + Cu(s)",
    components=components_2,
)

# NOTE: print analysis
print(
    f"[bold underline]Reaction Analysis for: {reaction_2.name}[/bold underline]")
print(f"Reaction: {reaction_2.reaction}")
print(f"Component IDs: {reaction_2.component_ids}")
print(f"Reaction Mode Symbol: {reaction_2.reaction_mode_symbol}")
print(
    f"reaction type: {reaction_2.reaction_type}"
)
print(
    f"Symbolic Unbalanced Reaction: {reaction_2.symbolic_unbalanced_reaction}")
print(f"Symbolic Reaction: {reaction_2.symbolic_reaction}")
print(f"Reactants: {reaction_2.reactants}")
print(f"Products: {reaction_2.products}")
print(f"Reaction Coefficients: {reaction_2.reaction_coefficients}")
print(f"Reaction Stoichiometry: {reaction_2.reaction_stoichiometry}")
print(
    f"Reaction Stoichiometry Matrix: {reaction_2.reaction_stoichiometry_matrix}")
print(
    f"Reaction Stoichiometry Source: {reaction_2.reaction_stoichiometry_source}")
print(f"State Counts: {reaction_2.state_count}")
print(f"Reaction Phase: {reaction_2.reaction_phase}")
print(f"Reaction State: {reaction_2.reaction_state}")
print(f"Carbon Count: {reaction_2.carbon_count}")
print(f"Reactants Names: {reaction_2.reactants_names}")
print(f"Products Names: {reaction_2.products_names}")
print(f"All Components: {reaction_2.all_components}")
print(f"Available components: {reaction_2.available_components}")
print(f"Component Checker: {reaction_2.component_checker}")
print(f"Mapped Components: {reaction_2.map_components}")

# NOTE: electron transfer reaction (electron requires an explicit state to
# parse; Fe{3+}/Fe{2+} share the "aq" state so no component list can uniquely
# map them either)
reaction_3 = Reaction(
    name="Iron(III) Reduction by Electron Transfer",
    reaction="Fe{3+}(aq) + e{-}(aq) => Fe{2+}(aq)",
    components=None,
)

# NOTE: print analysis
print(
    f"[bold underline]Reaction Analysis for: {reaction_3.name}[/bold underline]")
print(f"Reaction: {reaction_3.reaction}")
print(f"Component IDs: {reaction_3.component_ids}")
print(f"Reaction Mode Symbol: {reaction_3.reaction_mode_symbol}")
print(
    f"reaction type: {reaction_3.reaction_type}"
)
print(
    f"Symbolic Unbalanced Reaction: {reaction_3.symbolic_unbalanced_reaction}")
print(f"Symbolic Reaction: {reaction_3.symbolic_reaction}")
print(f"Reactants: {reaction_3.reactants}")
print(f"Products: {reaction_3.products}")
print(f"Reaction Coefficients: {reaction_3.reaction_coefficients}")
print(f"Reaction Stoichiometry: {reaction_3.reaction_stoichiometry}")
print(
    f"Reaction Stoichiometry Matrix: {reaction_3.reaction_stoichiometry_matrix}")
print(
    f"Reaction Stoichiometry Source: {reaction_3.reaction_stoichiometry_source}")
print(f"State Counts: {reaction_3.state_count}")
print(f"Reaction Phase: {reaction_3.reaction_phase}")
print(f"Reaction State: {reaction_3.reaction_state}")
print(f"Carbon Count: {reaction_3.carbon_count}")
print(f"Reactants Names: {reaction_3.reactants_names}")
print(f"Products Names: {reaction_3.products_names}")
print(f"All Components: {reaction_3.all_components}")
print(f"Available components: {reaction_3.available_components}")
print(f"Component Checker: {reaction_3.component_checker}")
print(f"Mapped Components: {reaction_3.map_components}")
