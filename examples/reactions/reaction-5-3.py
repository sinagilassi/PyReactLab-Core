# import libs
from pyreactlab_core.models.reaction import Reaction
# ! print
from examples.utils.tools import print_reaction_analysis
from pythermodb_settings.models import Component

# NOTE: carbon-containing acid/base reaction
component_ch3cooh = Component(
    name="Acetic Acid",
    formula="CH3COOH",
    state="aq"
)
component_h_ion = Component(
    name="Hydrogen Ion",
    formula="H{+}",
    state="aq"
)
component_ch3coo_ion = Component(
    name="Acetate Ion",
    formula="CH3COO{-}",
    state="aq"
)

components_3 = [
    component_ch3cooh,
    component_h_ion,
    component_ch3coo_ion,
]

reaction_3 = Reaction(
    name="Acetic Acid Dissociation",
    reaction="CH3COOH(aq) <=> H{+}(aq) + CH3COO{-}(aq)",
    components=components_3,
)

# NOTE: print analysis
print_reaction_analysis(reaction_3)
