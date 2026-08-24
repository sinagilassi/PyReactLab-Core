# import libs
from pyreactlab_core.models.reaction import Reaction
# ! print
from examples.utils.tools import print_reaction_analysis
from pythermodb_settings.models import Component

# NOTE: electron transfer reaction (electron requires an explicit state to
# parse; Fe{3+}/Fe{2+} share the "aq" state so no component list can uniquely
# map them either)
reaction_3 = Reaction(
    name="Iron(III) Reduction by Electron Transfer",
    reaction="Fe{3+}(aq) + e{-}(aq) => Fe{2+}(aq)",
    components=None,
)

# NOTE: print analysis
print_reaction_analysis(reaction_3)
