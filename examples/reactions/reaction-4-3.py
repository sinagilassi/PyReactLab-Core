# import libs
from pyreactlab_core.models.reaction import Reaction
# ! print
from examples.utils.tools import print_reaction_analysis
from pythermodb_settings.models import Component

# NOTE: trivalent ion reaction
component_fe3_ion = Component(
    name="Iron(III) Ion",
    formula="Fe{3+}",
    state="aq"
)
component_oh_ion = Component(
    name="Hydroxide Ion",
    formula="OH{-}",
    state="aq"
)
component_fe_oh_3 = Component(
    name="Iron(III) Hydroxide",
    formula="Fe(OH)3",
    state="s"
)

components_3 = [
    component_fe3_ion,
    component_oh_ion,
    component_fe_oh_3,
]

reaction_3 = Reaction(
    name="Iron(III) Hydroxide Precipitation",
    reaction="Fe{3+}(aq) + 3OH{-}(aq) => Fe(OH)3(s)",
    components=components_3,
)

# NOTE: print analysis
print_reaction_analysis(reaction_3)
