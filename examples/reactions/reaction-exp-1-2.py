# import libs
from pyreactlab_core.models.reaction import Reaction
from pythermodb_settings.models import Component
# ! print
from examples.utils.tools import print_reaction_analysis

# NOTE: define components
component_ch4 = Component(
    name="Methane",
    formula="CH4",
    state="l"
)
component_o2 = Component(
    name="Oxygen",
    formula="O2",
    state="g"
)

component_co2 = Component(
    name="Carbon Dioxide",
    formula="CO2",
    state="g"
)

component_h2o = Component(
    name="Water",
    formula="H2O",
    state="g"
)

components = [
    component_ch4,
    component_o2,
    component_co2,
    component_h2o
]

# components = []
# components = None

# NOTE: define components
reaction_1 = Reaction(
    name="Combustion of Methane",
    reaction="CH4(l) + 2O2(g) => CO2(g) + 2H2O(g)",
    components=components,
)

# NOTE: print analysis
print_reaction_analysis(reaction_1)
