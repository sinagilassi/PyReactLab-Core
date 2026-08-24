# import libs
from pyreactlab_core.models.reaction import Reaction
from pythermodb_settings.models import Component
from rich import print

# NOTE: precipitation reaction
component_ag_ion = Component(
    name="Silver Ion",
    formula="Ag",
    state="aq"
)
component_cl_ion = Component(
    name="Chloride Ion",
    formula="Cl",
    state="aq"
)
component_agcl = Component(
    name="Silver Chloride",
    formula="AgCl",
    state="s"
)

components_1 = [
    component_ag_ion,
    component_cl_ion,
    component_agcl,
]

reaction_1 = Reaction(
    name="Silver Chloride Precipitation",
    reaction="Ag{+}(aq) + Cl{-}(aq) => AgCl(s)",
    components=components_1,
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

# NOTE: divalent ions reaction
component_ca_ion = Component(
    name="Calcium Ion",
    formula="Ca",
    state="aq"
)
component_co3_ion = Component(
    name="Carbonate Ion",
    formula="CO3",
    state="aq"
)
component_caco3 = Component(
    name="Calcium Carbonate",
    formula="CaCO3",
    state="s"
)

components_2 = [
    component_ca_ion,
    component_co3_ion,
    component_caco3,
]

reaction_2 = Reaction(
    name="Calcium Carbonate Precipitation",
    reaction="Ca{2+}(aq) + CO3{2-}(aq) => CaCO3(s)",
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

# NOTE: trivalent ion reaction
component_fe3_ion = Component(
    name="Iron(III) Ion",
    formula="Fe",
    state="aq"
)
component_oh_ion = Component(
    name="Hydroxide Ion",
    formula="OH",
    state="aq"
)
component_fe_oh_3 = Component(
    name="Iron(III) Hydroxide",
    formula="Fe(OH)3",
    state="s"
)

components_3 = [
    component_fe3_ion,
    component_oh_ion,
    component_fe_oh_3,
]

reaction_3 = Reaction(
    name="Iron(III) Hydroxide Precipitation",
    reaction="Fe{3+}(aq) + 3OH{-}(aq) => Fe(OH)3(s)",
    components=components_3,
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
