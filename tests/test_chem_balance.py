from pyreactlab_core.docs.chem_balance import (
    _split_equation,
    balance,
    parse_species,
)


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
        "CO2 + 3 H2 -> CH3OH + H2O"
    )
    assert balance("CO2(g) + H2(g) => 3CH3OH(g) + H2O(g)") == (
        "CO2(g) + 3 H2(g) -> CH3OH(g) + H2O(g)"
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
