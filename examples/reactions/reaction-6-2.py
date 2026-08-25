# import libs
from pyreactlab_core.models.reaction import Reaction
# ! print
from examples.utils.tools import print_reaction_analysis
from pythermodb_settings.models import Component

# NOTE: redox reaction (distinct states keep Zn/Cu identities unique)
component_zn_s = Component(
    name="Zinc",
    formula="Zn",
    state="s"
)
component_cu_ion = Component(
    name="Copper(II) Ion",
    formula="Cu{2+}",
    state="aq"
)
component_zn_ion = Component(
    name="Zinc Ion",
    formula="Zn{2+}",
    state="aq"
)
component_cu_s = Component(
    name="Copper",
    formula="Cu",
    state="s"
)

components_2 = [
    component_zn_s,
    component_cu_ion,
    component_zn_ion,
    component_cu_s,
]

reaction_2 = Reaction(
    name="Zinc-Copper Redox Reaction",
    reaction="Zn(s) + Cu{2+}(aq) => Zn{2+}(aq) + Cu(s)",
    components=components_2,
)

# NOTE: print analysis
print_reaction_analysis(reaction_2)
