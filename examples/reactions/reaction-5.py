# import libs
from pyreactlab_core.models.reaction import Reaction
# ! print
from examples.utils.tools import print_reaction_analysis
from pythermodb_settings.models import Component

# NOTE: polyatomic ions reaction
component_ca_ion = Component(
    name="Calcium Ion",
    formula="Ca",
    state="aq"
)
component_po4_ion = Component(
    name="Phosphate Ion",
    formula="PO4",
    state="aq"
)
component_ca3_po4_2 = Component(
    name="Calcium Phosphate",
    formula="Ca3(PO4)2",
    state="s"
)

components_1 = [
    component_ca_ion,
    component_po4_ion,
    component_ca3_po4_2,
]

reaction_1 = Reaction(
    name="Calcium Phosphate Precipitation",
    reaction="3Ca{2+}(aq) + 2PO4{3-}(aq) => Ca3(PO4)2(s)",
    components=components_1,
)

# NOTE: print analysis
print_reaction_analysis(reaction_1)
