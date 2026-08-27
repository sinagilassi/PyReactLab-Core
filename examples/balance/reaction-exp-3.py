# import libs
from pyreactlab_core.models.reaction import Reaction
from pyreactlab_core.docs import balance
from rich import print

reaction_1 = Reaction(
    name="Combustion of Methane",
    reaction="100CO2(g) + H2(g) => 3CH3OH(g) + H2O(g)"
)

reaction_2 = Reaction(
    name="Combustion of Ethanol",
    reaction="C2H5OH(l) + O2(g) <=> CO2(g) + H2O(g)"
)

# "KMnO4 + HCl = KCl + MnCl2 + H2O + Cl2"
reaction_3 = Reaction(
    name="Redox Reaction of Potassium Permanganate and Hydrochloric Acid",
    reaction="KMnO4(g) + HCl(l) = KCl(l) + MnCl2(s) + H2O(g) + Cl2(l)"
)

# "C6H5COOH + O2 = CO2 + H2O"
reaction_4 = Reaction(
    name="Combustion of Benzoic Acid",
    reaction="C6H5COOH(l) + O2(g) => CO2(g) + H2O(l)"
)

# "CuSO4*5H2O = CuSO4 + H2O"
reaction_5 = Reaction(
    name="Dehydration of Copper(II) Sulfate Pentahydrate",
    reaction="CuSO4*5H2O(s) <=> CuSO4(l) + H2O(l)"
)

# Radical abstraction example
reaction_6 = "H2(g) + 5Cl{*}(g) => HCl(g) + H{*}(g)"

# Radical ion formation example
reaction_7 = "O2(g) + e- => 2O2{*-}(g)"

# Radical ion protonation example
reaction_8 = "O2{*-}(aq) + H{+}(aq) => HO2{*}(aq)"

# Zwitterion intramolecular proton transfer example
reaction_9 = "3NH3{+}-CH2-COO{-}(aq) => NH2-CH2-COOH(aq)"

# Zwitterion acid-base example with mixed charge centers
reaction_10 = (
    "NH3{+}-CH2-NH3{+}-COO{-}(aq) + OH{-}(aq) "
    "=> NH2-CH2-NH3{+}-COO{-}(aq) + 2H2O(l)"
)

tests = [
    reaction_1,
    reaction_2,
    reaction_3,
    reaction_4,
    reaction_5,
    reaction_6,
    reaction_7,
    reaction_8,
    reaction_9,
    reaction_10,
]

for eq in tests:
    reaction = eq.reaction if isinstance(eq, Reaction) else eq
    print("IN :", reaction)
    print("ALG:", balance(eq, "algebraic"))
    print("HAL:", balance(eq, "half", medium="auto"))
    print("OX :", balance(eq, "oxidation"))
    print("-" * 60)
