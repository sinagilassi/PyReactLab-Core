# import libs
from pyreactlab_core.docs import balance
from rich import print

tests = [
    "Fe + Cl2 = FeCl3",
    "KMnO4 + HCl = KCl + MnCl2 + H2O + Cl2",
    "K4Fe(CN)6 + H2SO4 + H2O = K2SO4 + FeSO4 + (NH4)2SO4 + CO",
    "C6H5COOH + O2 = CO2 + H2O",
    "CuSO4*5H2O = CuSO4 + H2O",
    # ionic-ish (you can use { } for charge like in your screenshot):
    "Cr2O7{2-} + H{+} + e- = Cr{3+} + H2O",
    "S{2-} + I2 = I{-} + S",
    "MnO4{-} + Fe{2+} + H{+} -> Mn{2+} + Fe{3+} + H2O",
    "Cu + HNO3 -> Cu(NO3)2 + NO2 + H2O",
    "Pb + PbO2 + H2SO4 -> PbSO4 + H2O",
    "FeS2 + O2 + H2O -> Fe(OH)3 + H2SO4",
    "MnO4{-} + C2O4{2-} + H{+} -> Mn{2+} + CO2 + H2O",
    "Au + CN{-} + O2 + H2O -> [Au(CN)2]{-} + OH{-}"
]


for eq in tests:
    print("IN :", eq)
    print("ALG:", balance(eq, "algebraic"))
    print("HAL:", balance(eq, "half", medium="auto"))
    print("OX :", balance(eq, "oxidation"))
    print("-" * 60)
