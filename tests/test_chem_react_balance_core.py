from pyreactlab_core.core.chem_react_balance import ChemReactBalance


def test_parse_ionic_charge_from_formula_notation():
    parser = ChemReactBalance()

    assert parser.parse_ionic_charge("Fe{3+}") == 3
    assert parser.parse_ionic_charge("SO4{2-}") == -2
    assert parser.parse_ionic_charge("OH{-}") == -1
    assert parser.parse_ionic_charge("Fe^3+") == 3
    assert parser.parse_ionic_charge("NH4+") == 1
    assert parser.parse_ionic_charge("Cl-") == -1
    assert parser.parse_ionic_charge("H2O") == 0
    assert parser.parse_ionic_charge("e{-}") == -1
    assert parser.parse_ionic_charge("CH3{*}") == 0
    assert parser.parse_ionic_charge("CH3{*+}") == 1
    assert parser.parse_ionic_charge("O2{*-}") == -1
    assert parser.parse_ionic_charge("C6H6{*2+}") == 2
    assert parser.parse_ionic_charge("NH3{+}-CH2-COO{-}") == 0
    assert parser.parse_ionic_charge("NH3{+}-CH2-NH3{+}-COO{-}") == 1


def test_parse_elemental_composition_ignores_supported_charge_notation():
    parser = ChemReactBalance()

    assert parser.parse_elemental_composition("SO4{2-}") == {
        "S": 1,
        "O": 4,
    }
    assert parser.parse_elemental_composition("Fe^3+") == {"Fe": 1}
    assert parser.parse_elemental_composition("NH4+") == {
        "N": 1,
        "H": 4,
    }
    assert parser.parse_elemental_composition("e{-}") == {}
    assert parser.parse_elemental_composition("CH3{*}") == {
        "C": 1,
        "H": 3,
    }
    assert parser.parse_elemental_composition("O2{*-}") == {
        "O": 2,
    }
    assert parser.parse_elemental_composition("NH3{+}-CH2-COO{-}") == {
        "N": 1,
        "H": 5,
        "C": 2,
        "O": 2,
    }
    assert parser.parse_elemental_composition("NH3{+}-CH2-NH3{+}-COO{-}") == {
        "N": 2,
        "H": 8,
        "C": 2,
        "O": 2,
    }
