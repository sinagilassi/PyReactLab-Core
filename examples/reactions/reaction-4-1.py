# import libs
from pyreactlab_core.models.reaction import Reaction
# ! print
from examples.utils.tools import print_reaction_analysis
from pythermodb_settings.models import Component

# NOTE: precipitation reaction
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
print_reaction_analysis(reaction_1)
