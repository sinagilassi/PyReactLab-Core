# import libs
from pyreactlab_core.models.reaction import Reaction
from pythermodb_settings.models import Component
from rich import print
from pyreactlab_core import build_rxns_stoichiometry

# NOTE: define ionic components
component_ag_ion = Component(
    name="Silver Ion",
    formula="Ag{+}",
    state="aq"
)

component_cl_ion = Component(
    name="Chloride Ion",
    formula="Cl{-}",
    state="aq"
)

component_agcl = Component(
    name="Silver Chloride",
    formula="AgCl",
    state="s"
)

component_ca_ion = Component(
    name="Calcium Ion",
    formula="Ca{2+}",
    state="aq"
)

component_co3_ion = Component(
    name="Carbonate Ion",
    formula="CO3{2-}",
    state="aq"
)

component_caco3 = Component(
    name="Calcium Carbonate",
    formula="CaCO3",
    state="s"
)

components = [
    component_ag_ion,
    component_cl_ion,
    component_agcl,
    component_ca_ion,
    component_co3_ion,
    component_caco3
]

# NOTE: define ionic reaction string
reaction_1 = "Ag{+}(aq) + Cl{-}(aq) => AgCl(s)"
name_1 = "Silver Chloride Precipitation"
components_1 = [component_ag_ion, component_cl_ion, component_agcl]

# second ionic reaction
reaction_2 = "Ca{2+}(aq) + CO3{2-}(aq) => CaCO3(s)"
name_2 = "Calcium Carbonate Precipitation"
components_2 = [component_ca_ion, component_co3_ion, component_caco3]

# NOTE: define reactions
reaction_1 = Reaction(
    name=name_1,
    reaction=reaction_1,
    components=components_1
)

reaction_2 = Reaction(
    name=name_2,
    reaction=reaction_2,
    components=components_2
)

# SECTION: print analysis
# Reaction 1
print(
    f"[bold underline]Reaction Analysis for: {reaction_1.name}[/bold underline]")
print(f"Reaction: {reaction_1.reaction}")
print(f"Component IDs: {reaction_1.component_ids}")
print(f"Reaction Coefficients: {reaction_1.reaction_coefficients}")
print(f"Reaction Stoichiometry: {reaction_1.reaction_stoichiometry}")
print(
    f"Reaction Stoichiometry Matrix: {reaction_1.reaction_stoichiometry_matrix}")
print(
    f"Reaction Stoichiometry Source: {reaction_1.reaction_stoichiometry_source}")
print(f"Charge Count: {reaction_1.charge_count}")
print(f"Net Charge: {reaction_1.net_charge}")

# Reaction 2
print(
    f"\n[bold underline]Reaction Analysis for: {reaction_2.name}[/bold underline]")
print(f"Reaction: {reaction_2.reaction}")
print(f"Component IDs: {reaction_2.component_ids}")
print(f"Reaction Coefficients: {reaction_2.reaction_coefficients}")
print(f"Reaction Stoichiometry: {reaction_2.reaction_stoichiometry}")
print(
    f"Reaction Stoichiometry Matrix: {reaction_2.reaction_stoichiometry_matrix}")
print(
    f"Reaction Stoichiometry Source: {reaction_2.reaction_stoichiometry_source}")
print(f"Charge Count: {reaction_2.charge_count}")
print(f"Net Charge: {reaction_2.net_charge}")

# Components
print("\n[bold underline]Components List[/bold underline]")
print(components)


# SECTION: Get reaction stoichiometry matrix
stoichiometry_result = build_rxns_stoichiometry(
    reactions=[reaction_1, reaction_2],
    components=components,
    component_key="Name-Formula"
)
print(stoichiometry_result)
