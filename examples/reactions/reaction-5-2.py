# import libs
from pyreactlab_core.models.reaction import Reaction
# ! print
from examples.utils.tools import print_reaction_analysis
from pythermodb_settings.models import Component

# NOTE: carbon subscript reaction
component_ca_ion = Component(
    name="Calcium Ion",
    formula="Ca{2+}",
    state="aq"
)
component_c2o4_ion = Component(
    name="Oxalate Ion",
    formula="C2O4{2-}",
    state="aq"
)
component_caC2o4 = Component(
    name="Calcium Oxalate",
    formula="CaC2O4",
    state="s"
)

components_2 = [
    component_ca_ion,
    component_c2o4_ion,
    component_caC2o4,
]

reaction_2 = Reaction(
    name="Calcium Oxalate Precipitation",
    reaction="Ca{2+}(aq) + C2O4{2-}(aq) => CaC2O4(s)",
    components=components_2,
)

# NOTE: print analysis
print_reaction_analysis(reaction_2)
