from pyreactlab_core.configs.constants import ALL_REACTION_MODE_SYMBOLS
from pyreactlab_core.docs.chem_balance import (
    _split_equation,
    balance,
    parse_species,
)
from pyreactlab_core.models.reaction import Reaction


def test_trailing_plus_ion_is_not_split_as_separator():
    assert _split_equation("NH4+ + OH{-} -> NH3 + H2O") == (
        ["NH4+", "OH{-}"],
        ["NH3", "H2O"],
    )

    ammonium = parse_species("NH4+")
    assert ammonium.formula == "NH4"
    assert ammonium.charge == 1
    assert ammonium.atoms == {"N": 1, "H": 4}

    assert balance("NH4+ + OH{-} -> NH3 + H2O") == (
        "NH4+ + OH{-} -> NH3 + H2O"
    )


def test_leading_coefficients_are_stripped_from_raw_string_species():
    assert balance("CO2 + H2 => 3CH3OH + H2O") == (
        "CO2 + 3 H2 => CH3OH + H2O"
    )
    assert balance("CO2(g) + H2(g) => 3CH3OH(g) + H2O(g)") == (
        "CO2(g) + 3 H2(g) => CH3OH(g) + H2O(g)"
    )


def test_balance_reaction_object_preserves_states():
    reaction = Reaction(
        name="Hydrogen Combustion",
        reaction="H2(g) + O2(g) => H2O(g)",
    )

    assert balance(reaction) == "2 H2(g) + O2(g) => 2 H2O(g)"


def test_reaction_model_does_not_auto_balance_by_default():
    reaction = Reaction(
        name="Hydrogen Combustion",
        reaction="H2(g) + O2(g) => H2O(g)",
    )

    assert reaction.reaction == "H2(g) + O2(g) => H2O(g)"
    assert reaction.is_balanced is False


def test_reaction_model_can_auto_balance_before_analysis():
    reaction = Reaction(
        name="Hydrogen Combustion",
        reaction="H2(g) + O2(g) => H2O(g)",
        balance_reaction=True,
    )

    assert reaction.reaction == "2 H2(g) + O2(g) => 2 H2O(g)"
    assert reaction.reaction_stoichiometry == {
        "H2-g": -2.0,
        "O2-g": -1.0,
        "H2O-g": 2.0,
    }
    assert reaction.is_balanced is True


def test_balance_preserves_original_reaction_operator():
    for operator in ALL_REACTION_MODE_SYMBOLS:
        assert balance(f"H2 + O2 {operator} H2O") == (
            f"2 H2 + O2 {operator} 2 H2O"
        )


def test_underdetermined_balance_keeps_all_input_species():
    assert balance("C + O2 -> CO + CO2") == "3 C + 2 O2 -> 2 CO + CO2"


def test_method_aliases_are_accepted():
    expected = "2 H2 + O2 -> 2 H2O"
    for method in (
        "matrix",
        "half-reaction",
        "ion",
        "ion-electron",
        "oxidation-number",
        "ox",
    ):
        assert balance("H2 + O2 -> H2O", method) == expected


def test_half_reaction_prefers_complete_original_equation_when_possible():
    assert balance(
        "MnO4{-} + Fe{2+} + H{+} -> Mn{2+} + Fe{3+} + H2O",
        "half",
    ) == "MnO4{-} + 5 Fe{2+} + 8 H{+} -> Mn{2+} + 5 Fe{3+} + 4 H2O"


def test_balance_parser_supports_radicals_and_zwitterions():
    radical_anion = parse_species("O2{*-}(g)")
    assert radical_anion.formula == "O2(g)"
    assert radical_anion.charge == -1
    assert radical_anion.atoms == {"O": 2}

    zwitterion = parse_species("NH3{+}-CH2-COO{-}(aq)")
    assert zwitterion.formula == "NH3-CH2-COO(aq)"
    assert zwitterion.charge == 0
    assert zwitterion.atoms == {
        "N": 1,
        "H": 5,
        "C": 2,
        "O": 2,
    }

    assert balance("H2(g) + Cl{*}(g) => HCl(g) + H{*}(g)") == (
        "H2(g) + Cl{*}(g) => HCl(g) + H{*}(g)"
    )
    assert balance("O2(g) + e- => O2{*-}(g)") == (
        "O2(g) + e- => O2{*-}(g)"
    )
    assert balance("O2{*-}(aq) + H{+}(aq) => HO2{*}(aq)") == (
        "O2{*-}(aq) + H{+}(aq) => HO2{*}(aq)"
    )
    assert balance("NH3{+}-CH2-COO{-}(aq) => NH2-CH2-COOH(aq)") == (
        "NH3{+}-CH2-COO{-}(aq) => NH2-CH2-COOH(aq)"
    )
