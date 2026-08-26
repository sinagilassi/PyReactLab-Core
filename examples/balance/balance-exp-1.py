# import libs
from rich import print
from pyreactlab_core.core import parse_elemental_composition, parse_ionic_charge

# NOTE: Examples of parsing elemental composition

compounds = [
    "Fe(OH)3",
    "Ca3(PO4)2",
    "CuSO4*5H2O",
    "SO4{2-}",
    "e{-}",
    "H2O",
    "C6H12O6",
]

# SECTION: Parse elemental composition for each compound
for compound in compounds:
    res_ = parse_elemental_composition(compound)
    print("Compound:")
    print(f"  {compound}")
    print("Elemental Composition:")
    print(f"  {res_}")
    print("-" * 40)

# SECTION: Parse ionic charge for each compound
for compound in compounds:
    res_ = parse_ionic_charge(compound)
    print("Compound:")
    print(f"  {compound}")
    print("Ionic Charge:")
    print(f"  {res_}")
    print("-" * 40)
