import pytest
from pythermodb_settings.models import Component

from pyreactlab_core.models import Reaction, ReactionNetwork


def make_network(name: str, reactions: list[tuple[str, str]]) -> ReactionNetwork:
    return ReactionNetwork(
        name=name,
        reactions=[
            Reaction(name=reaction_name, reaction=reaction)
            for reaction_name, reaction in reactions
        ],
    )


def test_sequential_network_species_roles_and_rank():
    network = make_network(
        "sequential",
        [
            ("R1", "H2(g) => 2H(g)"),
            ("R2", "2H(g) + O(g) => H2O(g)"),
        ],
    )

    assert network.source_species == ["H2-g", "O-g"]
    assert network.intermediate_species == ["H-g"]
    assert network.sink_species == ["H2O-g"]
    assert network.stoichiometric_rank == 2
    assert network.stoichiometric_matrix_shape == (4, 2)


def test_parallel_network_species_roles_are_deterministic():
    network = make_network(
        "parallel",
        [
            ("R1", "H2(g) => 2H(g)"),
            ("R2", "H2(g) => H2(l)"),
        ],
    )

    assert network.source_species == ["H2-g"]
    assert network.sink_species == ["H-g", "H2-l"]
    assert network.component_ids == ["H2-g", "H-g", "H2-l"]


def test_dependent_reaction_detection_is_incremental_and_deterministic():
    network = make_network(
        "dependent",
        [
            ("R1", "H2(g) => 2H(g)"),
            ("R2", "O2(g) => 2O(g)"),
            ("R3", "H2(g) + O2(g) => 2H(g) + 2O(g)"),
        ],
    )

    assert network.reaction_count == 3
    assert network.stoichiometric_rank == 2
    assert network.independent_reactions == ["R1", "R2"]
    assert network.dependent_reactions == ["R3"]
    assert network.independent_reaction_count == 2
    assert network.dependent_reaction_count == 1


def test_methanol_synthesis_network_structure():
    network = make_network(
        "methanol-synthesis",
        [
            ("R1", "CO2(g) + 3H2(g) <=> CH3OH(g) + H2O(g)"),
            ("R2", "CO2(g) + H2(g) <=> CO(g) + H2O(g)"),
            ("R3", "CO(g) + 2H2(g) <=> CH3OH(g)"),
        ],
    )

    assert network.component_ids == [
        "CO2-g",
        "H2-g",
        "CO-g",
        "CH3OH-g",
        "H2O-g",
    ]
    assert network.stoichiometric_matrix == [
        [-1.0, -1.0, 0.0],
        [-3.0, -1.0, -2.0],
        [0.0, 1.0, -1.0],
        [1.0, 0.0, 1.0],
        [1.0, 1.0, 0.0],
    ]
    assert network.stoichiometric_matrix_dict == {
        "CO2-g": [-1.0, -1.0, 0.0],
        "H2-g": [-3.0, -1.0, -2.0],
        "CO-g": [0.0, 1.0, -1.0],
        "CH3OH-g": [1.0, 0.0, 1.0],
        "H2O-g": [1.0, 1.0, 0.0],
    }
    assert network.stoichiometric_rank == 2
    assert network.dependent_reactions == ["R3"]
    assert network.reaction_mode_count == {
        "equilibrium": 0,
        "reversible": 3,
        "irreversible": 0,
    }
    assert network.phases == ["g"]
    assert network.components_by_phase == {
        "g": ["CO2-g", "H2-g", "CO-g", "CH3OH-g", "H2O-g"]
    }


def test_stoichiometric_matrix_can_be_reordered_by_components():
    component_co2 = Component(name="Carbon Dioxide", formula="CO2", state="g")
    component_h2 = Component(name="Hydrogen", formula="H2", state="g")
    component_co = Component(name="Carbon Monoxide", formula="CO", state="g")
    component_ch3oh = Component(name="Methanol", formula="CH3OH", state="g")
    component_h2o = Component(name="Water", formula="H2O", state="g")

    network = ReactionNetwork(
        name="methanol-synthesis",
        reactions=[
            Reaction(
                name="R1",
                reaction="CO2(g) + 3H2(g) <=> CH3OH(g) + H2O(g)",
                components=[
                    component_co2,
                    component_h2,
                    component_ch3oh,
                    component_h2o,
                ],
            ),
            Reaction(
                name="R2",
                reaction="CO2(g) + H2(g) <=> CO(g) + H2O(g)",
                components=[component_co2, component_h2, component_co, component_h2o],
            ),
            Reaction(
                name="R3",
                reaction="CO(g) + 2H2(g) <=> CH3OH(g)",
                components=[component_co, component_h2, component_ch3oh],
            ),
        ],
    )
    reordered_components = [
        component_h2o,
        component_ch3oh,
        component_co,
        component_h2,
        component_co2,
    ]

    assert network.stoichiometric_matrix_by_components(reordered_components) == [
        [1.0, 1.0, 0.0],
        [1.0, 0.0, 1.0],
        [0.0, 1.0, -1.0],
        [-3.0, -1.0, -2.0],
        [-1.0, -1.0, 0.0],
    ]
    assert network.stoichiometric_matrix_dict_by_components(reordered_components) == {
        "H2O-g": [1.0, 1.0, 0.0],
        "CH3OH-g": [1.0, 0.0, 1.0],
        "CO-g": [0.0, 1.0, -1.0],
        "H2-g": [-3.0, -1.0, -2.0],
        "CO2-g": [-1.0, -1.0, 0.0],
    }


def test_configured_components_define_network_component_order():
    component_co2 = Component(name="Carbon Dioxide", formula="CO2", state="g")
    component_h2 = Component(name="Hydrogen", formula="H2", state="g")
    component_co = Component(name="Carbon Monoxide", formula="CO", state="g")
    component_ch3oh = Component(name="Methanol", formula="CH3OH", state="g")
    component_h2o = Component(name="Water", formula="H2O", state="g")

    network = ReactionNetwork(
        name="methanol-synthesis",
        reactions=[
            Reaction(
                name="R1",
                reaction="CO2(g) + 3H2(g) <=> CH3OH(g) + H2O(g)",
                components=[
                    component_co2,
                    component_h2,
                    component_ch3oh,
                    component_h2o,
                ],
            ),
            Reaction(
                name="R2",
                reaction="CO2(g) + H2(g) <=> CO(g) + H2O(g)",
                components=[
                    component_co2,
                    component_h2,
                    component_co,
                    component_h2o,
                ],
            ),
            Reaction(
                name="R3",
                reaction="CO(g) + 2H2(g) <=> CH3OH(g)",
                components=[component_co, component_h2, component_ch3oh],
            ),
        ],
        components=[
            component_h2o,
            component_ch3oh,
            component_co,
            component_h2,
            component_co2,
        ],
    )

    assert network.component_ids == [
        "H2O-g",
        "CH3OH-g",
        "CO-g",
        "H2-g",
        "CO2-g",
    ]
    assert network.stoichiometric_matrix == [
        [1.0, 1.0, 0.0],
        [1.0, 0.0, 1.0],
        [0.0, 1.0, -1.0],
        [-3.0, -1.0, -2.0],
        [-1.0, -1.0, 0.0],
    ]
    assert network.stoichiometric_matrix_dict == {
        "H2O-g": [1.0, 1.0, 0.0],
        "CH3OH-g": [1.0, 0.0, 1.0],
        "CO-g": [0.0, 1.0, -1.0],
        "H2-g": [-3.0, -1.0, -2.0],
        "CO2-g": [-1.0, -1.0, 0.0],
    }


def test_configured_components_must_not_be_empty():
    reactions = [Reaction(name="R1", reaction="H2(g) => 2H(g)")]

    with pytest.raises(
        ValueError,
        match="ReactionNetwork components must not be empty when provided.",
    ):
        ReactionNetwork(name="empty-components", reactions=reactions, components=[])


def test_stoichiometric_matrix_by_components_defaults_to_network_order():
    network = make_network(
        "sequential",
        [
            ("R1", "H2(g) => 2H(g)"),
            ("R2", "2H(g) + O(g) => H2O(g)"),
        ],
    )

    assert network.stoichiometric_matrix_by_components() == network.stoichiometric_matrix
    assert (
        network.stoichiometric_matrix_dict_by_components()
        == network.stoichiometric_matrix_dict
    )


def test_stoichiometric_matrix_reorder_rejects_unknown_components():
    network = make_network("simple", [("R1", "H2(g) => 2H(g)")])
    component_o2 = Component(name="Oxygen", formula="O2", state="g")

    with pytest.raises(
        ValueError,
        match="Components are not part of this network: O2-g",
    ):
        network.stoichiometric_matrix_by_components([component_o2])


def test_ionic_network_preserves_charge_and_phase_identity():
    network = make_network(
        "neutralization",
        [("R1", "H{+}(aq) + OH{-}(aq) => H2O(l)")],
    )

    assert network.component_ids == ["H{+}-aq", "OH{-}-aq", "H2O-l"]
    assert network.charge_vector == [1.0, -1.0, 0.0]
    assert network.charge_balance_vector == [0.0]
    assert network.is_charge_balanced is True
    assert network.is_element_balanced is True
    assert network.stoichiometric_matrix == [[-1.0], [-1.0], [1.0]]


def test_multivalent_ions_remain_distinct_components():
    network = make_network(
        "iron-reduction",
        [("R1", "Fe{3+}(aq) + e{-}(aq) => Fe{2+}(aq)")],
    )

    assert network.component_ids == ["Fe{3+}-aq", "e{-}-aq", "Fe{2+}-aq"]
    assert network.charge_vector == [3.0, -1.0, 2.0]
    assert network.charge_balance_vector == [0.0]


def test_same_formula_different_phase_is_separate_component():
    network = make_network(
        "phase-change",
        [("R1", "H2O(l) => H2O(g)")],
    )

    assert network.component_ids == ["H2O-l", "H2O-g"]
    assert network.phases == ["l", "g"]
    assert network.multiphase_reactions == ["R1"]


def test_mixed_reaction_modes_are_aggregated():
    network = make_network(
        "modes",
        [
            ("R1", "H2(g) => 2H(g)"),
            ("R2", "2H(g) <=> H2(g)"),
            ("R3", "H2(g) = H2(l)"),
        ],
    )

    assert network.reaction_modes == {
        "R1": "irreversible",
        "R2": "reversible",
        "R3": "equilibrium",
    }
    assert network.equilibrium_reactions == ["R3"]
    assert network.reversible_reactions == ["R2"]
    assert network.irreversible_reactions == ["R1"]
    assert network.reaction_mode_count == {
        "equilibrium": 1,
        "reversible": 1,
        "irreversible": 1,
    }


def test_unbalanced_reaction_is_inspectable():
    network = make_network(
        "unbalanced",
        [("R1", "H2(g) => H(g)")],
    )

    assert network.is_balanced is False
    assert network.unbalanced_reactions == ["R1"]
    assert network.unbalanced_element_reactions == ["R1"]
    assert network.summary["is_balanced"] is False


def test_duplicate_reaction_names_are_rejected():
    reactions = [
        Reaction(name="R1", reaction="H2(g) => 2H(g)"),
        Reaction(name="R1", reaction="O2(g) => 2O(g)"),
    ]

    with pytest.raises(ValueError, match="Duplicate reaction name detected: R1"):
        ReactionNetwork(name="duplicate", reactions=reactions)
