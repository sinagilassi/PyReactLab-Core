import numpy as np

from pyreactlab_core.models import Reaction, ReactionNetwork


def make_network(name: str, reactions: list[tuple[str, str]]) -> ReactionNetwork:
    return ReactionNetwork(
        name=name,
        reactions=[
            Reaction(name=reaction_name, reaction=reaction)
            for reaction_name, reaction in reactions
        ],
    )


def assert_dependency_reconstructs(
    network: ReactionNetwork,
    dependent_reaction: str,
) -> None:
    matrix = np.array(network.stoichiometric_matrix, dtype=float)
    dependent_index = network.reaction_ids.index(dependent_reaction)
    dependent_vector = matrix[:, dependent_index]
    reconstructed = np.zeros_like(dependent_vector)

    for reaction_id, coefficient in network.reaction_dependencies[
        dependent_reaction
    ].items():
        reaction_index = network.reaction_ids.index(reaction_id)
        reconstructed += coefficient * matrix[:, reaction_index]

    assert np.allclose(reconstructed, dependent_vector)


def test_all_reactions_independent_have_no_dependencies():
    network = make_network(
        "all-independent",
        [
            ("R1", "H2(g) => 2H(g)"),
            ("R2", "O2(g) => 2O(g)"),
            ("R3", "N2(g) => 2N(g)"),
        ],
    )

    assert network.stoichiometric_rank == 3
    assert network.independent_reaction_indices == [0, 1, 2]
    assert network.dependent_reaction_indices == []
    assert network.independent_reactions == ["R1", "R2", "R3"]
    assert network.dependent_reactions == []
    assert network.reaction_dependencies == {}


def test_additive_dependent_reaction_coefficients():
    network = make_network(
        "additive",
        [
            ("R1", "H2(g) => 2H(g)"),
            ("R2", "O2(g) => 2O(g)"),
            ("R3", "H2(g) + O2(g) => 2H(g) + 2O(g)"),
        ],
    )

    assert network.independent_reactions == ["R1", "R2"]
    assert network.dependent_reactions == ["R3"]
    assert network.reaction_dependencies == {
        "R3": {
            "R1": 1.0,
            "R2": 1.0,
        }
    }
    assert_dependency_reconstructs(network, "R3")


def test_duplicate_reaction_dependency():
    network = make_network(
        "duplicate-stoichiometry",
        [
            ("R1", "H2(g) => 2H(g)"),
            ("R2", "H2(g) => 2H(g)"),
        ],
    )

    assert network.independent_reactions == ["R1"]
    assert network.dependent_reactions == ["R2"]
    assert network.reaction_dependencies == {"R2": {"R1": 1.0}}
    assert_dependency_reconstructs(network, "R2")


def test_scaled_reaction_dependency():
    network = make_network(
        "scaled",
        [
            ("R1", "H2(g) => 2H(g)"),
            ("R2", "2H2(g) => 4H(g)"),
        ],
    )

    assert network.independent_reactions == ["R1"]
    assert network.dependent_reactions == ["R2"]
    assert network.reaction_dependencies == {"R2": {"R1": 2.0}}
    assert_dependency_reconstructs(network, "R2")


def test_reverse_reaction_dependency():
    network = make_network(
        "reverse",
        [
            ("R1", "H2(g) => 2H(g)"),
            ("R2", "2H(g) => H2(g)"),
        ],
    )

    assert network.independent_reactions == ["R1"]
    assert network.dependent_reactions == ["R2"]
    assert network.reaction_dependencies == {"R2": {"R1": -1.0}}
    assert_dependency_reconstructs(network, "R2")


def test_multiple_dependent_reactions():
    network = make_network(
        "multiple-dependent",
        [
            ("R1", "H2(g) => 2H(g)"),
            ("R2", "O2(g) => 2O(g)"),
            ("R3", "H2(g) + O2(g) => 2H(g) + 2O(g)"),
            ("R4", "2H2(g) => 4H(g)"),
            ("R5", "2H(g) + O2(g) => H2(g) + 2O(g)"),
        ],
    )

    assert network.independent_reactions == ["R1", "R2"]
    assert network.dependent_reactions == ["R3", "R4", "R5"]
    assert network.reaction_dependencies == {
        "R3": {
            "R1": 1.0,
            "R2": 1.0,
        },
        "R4": {
            "R1": 2.0,
        },
        "R5": {
            "R1": -1.0,
            "R2": 1.0,
        },
    }
    for reaction_id in network.dependent_reactions:
        assert_dependency_reconstructs(network, reaction_id)


def test_methanol_synthesis_dependency_coefficients():
    network = make_network(
        "methanol-synthesis",
        [
            ("R1", "CO2(g) + 3H2(g) <=> CH3OH(g) + H2O(g)"),
            ("R2", "CO2(g) + H2(g) <=> CO(g) + H2O(g)"),
            ("R3", "CO(g) + 2H2(g) <=> CH3OH(g)"),
        ],
    )

    assert network.stoichiometric_rank == 2
    assert network.independent_reaction_indices == [0, 1]
    assert network.dependent_reaction_indices == [2]
    assert network.independent_reactions == ["R1", "R2"]
    assert network.dependent_reactions == ["R3"]
    assert network.reaction_dependencies == {
        "R3": {
            "R1": 1.0,
            "R2": -1.0,
        }
    }
    assert_dependency_reconstructs(network, "R3")


def test_deterministic_basis_follows_reaction_order():
    first_order = make_network(
        "first-order",
        [
            ("R1", "H2(g) => 2H(g)"),
            ("R2", "O2(g) => 2O(g)"),
            ("R3", "H2(g) + O2(g) => 2H(g) + 2O(g)"),
        ],
    )
    second_order = make_network(
        "second-order",
        [
            ("R3", "H2(g) + O2(g) => 2H(g) + 2O(g)"),
            ("R1", "H2(g) => 2H(g)"),
            ("R2", "O2(g) => 2O(g)"),
        ],
    )

    assert first_order.independent_reactions == ["R1", "R2"]
    assert second_order.independent_reactions == ["R3", "R1"]
    assert first_order.dependent_reactions == ["R3"]
    assert second_order.dependent_reactions == ["R2"]


def test_ionic_network_has_no_dependency_when_single_reaction():
    network = make_network(
        "ionic",
        [("R1", "H{+}(aq) + OH{-}(aq) => H2O(l)")],
    )

    assert network.independent_reactions == ["R1"]
    assert network.dependent_reactions == []
    assert network.reaction_dependencies == {}


def test_multiphase_duplicate_dependency_preserves_phase_identity():
    network = make_network(
        "multiphase-duplicate",
        [
            ("R1", "H2O(l) => H2O(g)"),
            ("R2", "2H2O(l) => 2H2O(g)"),
        ],
    )

    assert network.component_ids == ["H2O-l", "H2O-g"]
    assert network.independent_reactions == ["R1"]
    assert network.dependent_reactions == ["R2"]
    assert network.reaction_dependencies == {"R2": {"R1": 2.0}}
    assert_dependency_reconstructs(network, "R2")
