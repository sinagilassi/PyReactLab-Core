# import libs
from pyreactlab_core.models.reaction import Reaction
from pythermodb_settings.models import Component
# ! print
from examples.utils.tools import print_reaction_analysis

# NOTE: define components
component_co2 = Component(
    name="Carbon Dioxide",
    formula="CO2",
    state="g"
)

component_h2 = Component(
    name="Hydrogen",
    formula="H2",
    state="g"
)

component_ch3oh = Component(
    name="Methanol",
    formula="CH3OH",
    state="g"
)

component_h2o = Component(
    name="Water",
    formula="H2O",
    state="g"
)

component_co = Component(
    name="Carbon Monoxide",
    formula="CO",
    state="g"
)

components = [
    component_h2,
    component_co,
    component_ch3oh,
    component_h2o,
    component_co2,
]

# components = []
# components = None

# NOTE: define the reaction
# ! no balancing
reaction_1 = Reaction(
    name="Combustion of Methane",
    reaction="7CO2(g) + 3H2(g) ⇌ 2CH3OH(g) + H2O(g)",
    balance_reaction=False,
    components=components,
)

# NOTE: print analysis
print_reaction_analysis(reaction_1)

# ! balancing
reaction_2 = Reaction(
    name="Combustion of Methane",
    reaction="7CO2(g) + 3H2(g) ⇌ 2CH3OH(g) + H2O(g)",
    balance_reaction=True,
    components=components,
)

# NOTE: print analysis
print_reaction_analysis(reaction_2)
