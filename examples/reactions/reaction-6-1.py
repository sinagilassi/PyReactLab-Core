# import libs
from pyreactlab_core.models.reaction import Reaction
# ! print
from examples.utils.tools import print_reaction_analysis
from pythermodb_settings.models import Component

# NOTE: C-starting element, but NOT carbon (Fe/Ce both change oxidation state,
# so the plain-formula identity used for charged species collides and no
# component list can uniquely map them)
reaction_1 = Reaction(
    name="Cerium-Iron Redox Reaction",
    reaction="Fe{2+}(aq) + Ce{4+}(aq) => Fe{3+}(aq) + Ce{3+}(aq)",
    components=None,
)

# NOTE: print analysis
print_reaction_analysis(reaction_1)
