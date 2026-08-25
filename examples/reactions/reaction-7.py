from rich import print
from rich.table import Table

from pyreactlab_core.models.reaction import Reaction


reaction_cases = [
    {
        "name": "Hydrogen Combustion",
        "reaction": "2H2(g) + O2(g) => 2H2O(g)",
        "expected": (True, True, True),
    },
    {
        "name": "Unbalanced Hydrogen Combustion",
        "reaction": "H2(g) + O2(g) => H2O(g)",
        "expected": (False, True, False),
    },
    {
        "name": "Iron(III) Hydroxide Precipitation",
        "reaction": "Fe{3+}(aq) + 3OH{-}(aq) => Fe(OH)3(s)",
        "expected": (True, True, True),
    },
    {
        "name": "Neutralization",
        "reaction": "H{+}(aq) + OH{-}(aq) => H2O(l)",
        "expected": (True, True, True),
    },
    {
        "name": "Sulfate Charge Loss",
        "reaction": "SO4{2-}(aq) => SO4(aq)",
        "expected": (True, False, False),
    },
    {
        "name": "Calcium Carbonate Decomposition",
        "reaction": "CaCO3(s) => CaO(s) + CO2(g)",
        "expected": (True, True, True),
    },
    {
        "name": "Potassium Chlorate Decomposition",
        "reaction": "2KClO3(s) => 2KCl(s) + 3O2(g)",
        "expected": (True, True, True),
    },
    {
        "name": "Copper Sulfate Pentahydrate Dehydration",
        "reaction": "CuSO4*5H2O(s) => CuSO4(s) + 5H2O(l)",
        "expected": (True, True, True),
    },
    {
        "name": "Calcium Phosphate Decomposition",
        "reaction": "Ca3(PO4)2(s) => 3CaO(s) + P2O5(s)",
        "expected": (True, True, True),
    },
]


def status(value: bool) -> str:
    return "[green]PASS[/green]" if value else "[red]FAIL[/red]"


table = Table(title="Element / Charge / Overall Balance")
table.add_column("Reaction")
table.add_column("Element", justify="center")
table.add_column("Charge", justify="center")
table.add_column("Overall", justify="center")

for case in reaction_cases:
    reaction = Reaction(
        name=case["name"],
        reaction=case["reaction"],
        components=None,
    )

    actual = (
        reaction.is_element_balanced,
        reaction.is_charge_balanced,
        reaction.is_balanced,
    )

    assert actual == case["expected"], (
        f"{case['reaction']} expected {case['expected']} but got {actual}"
    )

    table.add_row(
        reaction.reaction,
        status(reaction.is_element_balanced),
        status(reaction.is_charge_balanced),
        status(reaction.is_balanced),
    )

print(table)
