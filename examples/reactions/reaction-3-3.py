# import libs
from pyreactlab_core.models.reaction import Reaction
# ! print
from examples.utils.tools import print_reaction_analysis
from pythermodb_settings.models import Component

# NOTE: aqueous neutralization reaction
component_h_ion = Component(
    name="Hydrogen Ion",
    formula="H{+}",
    state="aq"
)
component_oh_ion = Component(
    name="Hydroxide Ion",
    formula="OH{-}",
    state="aq"
)
component_h2o_l = Component(
    name="Water",
    formula="H2O",
    state="l"
)

components_3 = [
    component_h_ion,
    component_oh_ion,
    component_h2o_l,
]

reaction_3 = Reaction(
    name="Acid-Base Neutralization",
    reaction="H{+}(aq) + OH{-}(aq) => H2O(l)",
    components=components_3,
)

# NOTE: print analysis
print_reaction_analysis(reaction_3)
