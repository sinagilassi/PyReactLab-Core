# import libs
from pyreactlab_core.models.reaction import Reaction
# ! print
from examples.utils.tools import print_reaction_analysis
from pythermodb_settings.models import Component

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
print_reaction_analysis(reaction_2)
