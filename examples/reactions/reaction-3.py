# import libs
from pyreactlab_core.models.reaction import Reaction
# ! print
from examples.utils.tools import print_reaction_analysis
from pythermodb_settings.models import Component

# NOTE: gas phase reaction
component_h2 = Component(
    name="Hydrogen",
    formula="H2",
    state="g"
)
component_o2 = Component(
    name="Oxygen",
    formula="O2",
    state="g"
)
component_h2o_g = Component(
    name="Water",
    formula="H2O",
    state="g"
)

components_1 = [
    component_h2,
    component_o2,
    component_h2o_g,
]

reaction_1 = Reaction(
    name="Combustion of Hydrogen",
    reaction="2H2(g) + O2(g) => 2H2O(g)",
    components=components_1,
)

# NOTE: print analysis
print_reaction_analysis(reaction_1)
