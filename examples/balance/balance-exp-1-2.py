# import libs
from rich import print
from pyreactlab_core.core import parse_elemental_composition, parse_ionic_charge

# NOTE: Examples of parsing elemental composition

CASES = [
    {
        "name": "Water",
        "formula": "H2O",
        "state": "l",
        "charge": 0,
        "species_type": ["neutral"],
    },
    {
        "name": "Dioxygen",
        "formula": "O2",
        "state": "g",
        "charge": 0,
        "species_type": ["neutral"],
    },
    {
        "name": "Carbon dioxide",
        "formula": "CO2",
        "state": "g",
        "charge": 0,
        "species_type": ["neutral"],
    },
    {
        "name": "Methanol",
        "formula": "CH4O",
        "state": "l",
        "charge": 0,
        "species_type": ["neutral"],
    },
    {
        "name": "Benzene",
        "formula": "C6H6",
        "state": "l",
        "charge": 0,
        "species_type": ["neutral"],
    },
    {
        "name": "Iron(III)",
        "formula": "Fe{3+}",
        "state": "s",
        "charge": 3,
        "species_type": ["cation"],
    },
    {
        "name": "Copper(II)",
        "formula": "Cu{2+}",
        "state": "aq",
        "charge": 2,
        "species_type": ["cation"],
    },
    {
        "name": "Calcium",
        "formula": "Ca{2+}",
        "state": "aq",
        "charge": 2,
        "species_type": ["cation"],
    },
    {
        "name": "Aluminum",
        "formula": "Al{3+}",
        "state": "aq",
        "charge": 3,
        "species_type": ["cation"],
    },
    {
        "name": "Cerium(IV)",
        "formula": "Ce{4+}",
        "state": "aq",
        "charge": 4,
        "species_type": ["cation"],
    },
    {
        "name": "Ammonium",
        "formula": "NH4{+}",
        "state": "aq",
        "charge": 1,
        "species_type": ["cation"],
    },
    {
        "name": "Sulfate",
        "formula": "SO4{2-}",
        "state": "aq",
        "charge": -2,
        "species_type": ["anion"],
    },
    {
        "name": "Chloride",
        "formula": "Cl{-}",
        "state": "aq",
        "charge": -1,
        "species_type": ["anion"],
    },
    {
        "name": "Nitrate",
        "formula": "NO3{-}",
        "state": "aq",
        "charge": -1,
        "species_type": ["anion"],
    },
    {
        "name": "Carbonate",
        "formula": "CO3{2-}",
        "state": "aq",
        "charge": -2,
        "species_type": ["anion"],
    },
    {
        "name": "Phosphate",
        "formula": "PO4{3-}",
        "state": "aq",
        "charge": -3,
        "species_type": ["anion"],
    },
    {
        "name": "Ferricyanide",
        "formula": "Fe(CN)6{3-}",
        "state": "aq",
        "charge": -3,
        "species_type": ["anion"],
    },
    {
        "name": "Bromide",
        "formula": "Br{-}",
        "state": "aq",
        "charge": -1,
        "species_type": ["anion"],
    },
    {
        "name": "Methyl radical",
        "formula": "CH3{*}",
        "state": "g",
        "charge": 0,
        "species_type": ["neutral", "radical"],
    },
    {
        "name": "Hydroxyl radical",
        "formula": "HO{*}",
        "state": "g",
        "charge": 0,
        "species_type": ["neutral", "radical"],
    },
    {
        "name": "Ethyl radical",
        "formula": "C2H5{*}",
        "state": "g",
        "charge": 0,
        "species_type": ["neutral", "radical"],
    },
    {
        "name": "Nitric oxide radical",
        "formula": "NO{*}",
        "state": "g",
        "charge": 0,
        "species_type": ["neutral", "radical"],
    },
    {
        "name": "Methyl radical cation",
        "formula": "CH3{*+}",
        "state": "g",
        "charge": 1,
        "species_type": ["cation", "radical"],
    },
    {
        "name": "Superoxide radical anion",
        "formula": "O2{*-}",
        "state": "aq",
        "charge": -1,
        "species_type": ["anion", "radical"],
    },
    {
        "name": "Benzene radical cation",
        "formula": "C6H6{*+}",
        "state": "g",
        "charge": 1,
        "species_type": ["cation", "radical"],
    },
    {
        "name": "Naphthalene radical anion",
        "formula": "C10H8{*-}",
        "state": "g",
        "charge": -1,
        "species_type": ["anion", "radical"],
    },
    {
        "name": "Radical dication",
        "formula": "C6H6{*2+}",
        "state": "g",
        "charge": 2,
        "species_type": ["cation", "radical"],
    },
    {
        "name": "Radical trication",
        "formula": "Fe{*3+}",
        "state": "s",
        "charge": 3,
        "species_type": ["cation", "radical"],
    },
    {
        "name": "Radical dianion",
        "formula": "S{*2-}",
        "state": "s",
        "charge": -2,
        "species_type": ["anion", "radical"],
    },
    {
        "name": "Radical tetracation",
        "formula": "Ce{*4+}",
        "state": "s",
        "charge": 4,
        "species_type": ["cation", "radical"],
    },
    {
        "name": "Radical pentacation",
        "formula": "Ta{*5+}",
        "state": "s",
        "charge": 5,
        "species_type": ["cation", "radical"],
    },
    {
        "name": "Radical trianion",
        "formula": "P{*3-}",
        "state": "s",
        "charge": -3,
        "species_type": ["anion", "radical"],
    },
    {
        "name": "Radical tetraanion",
        "formula": "S{*4-}",
        "state": "s",
        "charge": -4,
        "species_type": ["anion", "radical"],
    },
    {
        "name": "Peroxide radical dianion",
        "formula": "O2{*2-}",
        "state": "aq",
        "charge": -2,
        "species_type": ["anion", "radical"],
    },
    {
        "name": "Metal radical dication",
        "formula": "Fe{*2+}",
        "state": "aq",
        "charge": 2,
        "species_type": ["cation", "radical"],
    },
    {
        "name": "Glycine zwitterion",
        "formula": "NH3{+}-CH2-COO{-}",
        "state": "s",
        "charge": 0,
        "species_type": ["zwitterion"],
    },
    {
        "name": "Alanine zwitterion",
        "formula": "NH3{+}-CH(CH3)-COO{-}",
        "state": "s",
        "charge": 0,
        "species_type": ["zwitterion"],
    },
    {
        "name": "Betaine zwitterion",
        "formula": "N(CH3)3{+}-CH2-COO{-}",
        "state": "s",
        "charge": 0,
        "species_type": ["zwitterion"],
    },
    {
        "name": "Aspartate dianion",
        "formula": "OOC{-}-CH2-CH(NH2)-COO{-}",
        "state": "aq",
        "charge": -2,
        "species_type": ["anion"],
    },
    {
        "name": "Net neutral mixed centers",
        "formula": "Fe{2+}-Cl{-}-Cl{-}",
        "state": "s",
        "charge": 0,
        "species_type": ["zwitterion"],
    },
    {
        "name": "Net dianion with mixed centers",
        "formula": "Fe{2+}-O{2-}-O{2-}",
        "state": "s",
        "charge": -2,
        "species_type": ["anion"],
    },
    {
        "name": "Net cation with mixed charge centers",
        "formula": "NH3{+}-CH2-NH3{+}-COO{-}",
        "state": "aq",
        "charge": 1,
        "species_type": ["cation"],
    },
    {
        "name": "Net dication with mixed centers",
        "formula": "Fe{3+}-Cl{-}",
        "state": "s",
        "charge": 2,
        "species_type": ["cation"],
    },
    {
        "name": "Radical zwitterion",
        "formula": "N{+}-O{-}-C{*}",
        "state": "s",
        "charge": 0,
        "species_type": ["zwitterion", "radical"],
    },
    {
        "name": "Radical net cation mixed centers",
        "formula": "N{+}-O{-}-C{*+}",
        "state": "s",
        "charge": 1,
        "species_type": ["cation", "radical"],
    },
    {
        "name": "Radical net dication mixed centers",
        "formula": "Fe{2+}-Cl{-}-C{*+}",
        "state": "s",
        "charge": 2,
        "species_type": ["cation", "radical"],
    },
    {
        "name": "Radical net trication mixed centers",
        "formula": "Fe{3+}-Cl{-}-C{*+}",
        "state": "s",
        "charge": 3,
        "species_type": ["cation", "radical"],
    },
    {
        "name": "Radical net dianion mixed centers",
        "formula": "Na{+}-O{2-}-O{*-}",
        "state": "s",
        "charge": -2,
        "species_type": ["anion", "radical"],
    },
    {
        "name": "Radical net trianion mixed centers",
        "formula": "Na{+}-P{3-}-O{*-}",
        "state": "s",
        "charge": -3,
        "species_type": ["anion", "radical"],
    },
]

# SECTION: Parse elemental composition for each compound
for compound in CASES:
    # compound
    compound = compound["formula"]
    res_ = parse_elemental_composition(compound)
    print("Compound:")
    print(f"  {compound}")
    print("Elemental Composition:")
    print(f"  {res_}")
    print("-" * 40)

# SECTION: Parse ionic charge for each compound
for compound in CASES:
    compound = compound["formula"]
    res_ = parse_ionic_charge(compound)
    print("Compound:")
    print(f"  {compound}")
    print("Ionic Charge:")
    print(f"  {res_}")
    print("-" * 40)
