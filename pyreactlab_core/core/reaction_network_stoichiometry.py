from __future__ import annotations

# import libs
import logging
from typing import Any

import numpy as np


# NOTE: set up logger
logger = logging.getLogger(__name__)

# NOTE: shared matrix type aliases
StoichiometricMatrix = list[list[float]]
ParticipationMatrix = list[list[int]]
ReactionDependencyMap = dict[str, dict[str, float]]

# NOTE: numerical tolerances for dependency reconstruction
DEFAULT_DEPENDENCY_ATOL = 1e-10
DEFAULT_DEPENDENCY_RTOL = 1e-9


# SECTION: stoichiometric matrix and rank analysis
def stoichiometric_matrix(
    reactions: list[Any],
    component_ids: list[str],
) -> StoichiometricMatrix:
    """
    Build the stoichiometric matrix using rows as components and columns as reactions.
    """
    # NOTE: reactant/product signs come from Reaction.reaction_stoichiometry.
    matrix: StoichiometricMatrix = []

    # SECTION: build matrix rows by component id
    for component_id in component_ids:
        row: list[float] = []
        for reaction in reactions:
            row.append(float(reaction.reaction_stoichiometry.get(component_id, 0.0)))
        matrix.append(row)

    return matrix


def stoichiometric_rank(matrix: StoichiometricMatrix) -> int:
    """
    Return the numeric rank of a stoichiometric matrix.
    """
    # ! empty matrices have no stoichiometric rank.
    if not matrix or not matrix[0]:
        return 0

    # NOTE: NumPy is already a package dependency.
    return int(np.linalg.matrix_rank(np.array(matrix, dtype=float)))


def clean_dependency_coefficient(
    value: float,
    tolerance: float = DEFAULT_DEPENDENCY_ATOL,
) -> float:
    """
    Remove insignificant numerical noise from dependency coefficients.
    """
    # NOTE: coefficients very close to zero are omitted by callers.
    if abs(value) < tolerance:
        return 0.0

    # NOTE: preserve simple integer stoichiometric relationships.
    nearest_integer = round(value)
    if abs(value - nearest_integer) < tolerance:
        return float(nearest_integer)

    return float(value)


def independent_reaction_indices(
    matrix: StoichiometricMatrix,
) -> list[int]:
    """
    Return deterministic independent reaction column indices.
    """
    # ! empty matrices cannot identify independent reactions.
    if not matrix or not matrix[0]:
        return []

    # SECTION: initialize incremental rank state
    array = np.array(matrix, dtype=float)
    selected_indexes: list[int] = []
    current_rank = 0

    # SECTION: scan reaction columns in user-provided order
    for index in range(array.shape[1]):
        candidate_indexes = selected_indexes + [index]
        candidate_rank = int(np.linalg.matrix_rank(array[:, candidate_indexes]))

        # NOTE: a reaction is independent when it increases matrix rank.
        if candidate_rank > current_rank:
            selected_indexes.append(index)
            current_rank = candidate_rank

    return selected_indexes


def dependent_reaction_indices(
    matrix: StoichiometricMatrix,
    independent_indices_value: list[int] | None = None,
) -> list[int]:
    """
    Return deterministic dependent reaction column indices.
    """
    # ! empty matrices cannot identify dependent reactions.
    if not matrix or not matrix[0]:
        return []

    # SECTION: determine independent indices when not supplied
    if independent_indices_value is None:
        independent_indices_value = independent_reaction_indices(matrix)

    independent_set = set(independent_indices_value)
    reaction_count = len(matrix[0])

    # NOTE: preserve original reaction column order.
    return [
        index
        for index in range(reaction_count)
        if index not in independent_set
    ]


def independent_stoichiometric_matrix(
    matrix: StoichiometricMatrix,
    independent_indices_value: list[int],
) -> StoichiometricMatrix:
    """
    Extract the selected independent reaction columns from a matrix.
    """
    # ! no independent columns means the basis matrix is empty.
    if not matrix or not independent_indices_value:
        return [[] for _ in matrix]

    array = np.array(matrix, dtype=float)
    return array[:, independent_indices_value].tolist()


def reaction_dependency_coefficients_by_index(
    matrix: StoichiometricMatrix,
    independent_indices_value: list[int],
    dependent_index: int,
    atol: float = DEFAULT_DEPENDENCY_ATOL,
    rtol: float = DEFAULT_DEPENDENCY_RTOL,
) -> dict[int, float]:
    """
    Solve one dependent reaction column against the independent reaction basis.
    """
    # ! dependency coefficients require a non-empty stoichiometric matrix.
    if not matrix or not matrix[0]:
        return {}

    array = np.array(matrix, dtype=float)
    dependent_vector = array[:, dependent_index]

    # ! a zero vector is dependent on any basis with all coefficients omitted.
    if np.allclose(dependent_vector, 0.0, atol=atol, rtol=rtol):
        return {}

    if not independent_indices_value:
        raise ValueError(
            "Cannot calculate dependency coefficients without an independent basis."
        )

    # SECTION: solve least-squares reconstruction against selected basis
    basis = array[:, independent_indices_value]
    coefficients, _, _, _ = np.linalg.lstsq(
        basis,
        dependent_vector,
        rcond=None,
    )

    reconstructed = basis @ coefficients

    # ! do not expose dependency coefficients unless reconstruction is valid.
    if not np.allclose(
        reconstructed,
        dependent_vector,
        atol=atol,
        rtol=rtol,
    ):
        raise ValueError(
            f"Dependent reaction column {dependent_index} could not be "
            "reconstructed from the independent basis."
        )

    # SECTION: clean coefficients and omit numerical zeros
    cleaned: dict[int, float] = {}
    for index, coefficient in zip(independent_indices_value, coefficients):
        cleaned_coefficient = clean_dependency_coefficient(
            float(coefficient),
            tolerance=atol,
        )
        if cleaned_coefficient != 0.0:
            cleaned[index] = cleaned_coefficient

    return cleaned


def reaction_dependencies(
    matrix: StoichiometricMatrix,
    ids: list[str],
    atol: float = DEFAULT_DEPENDENCY_ATOL,
    rtol: float = DEFAULT_DEPENDENCY_RTOL,
) -> ReactionDependencyMap:
    """
    Return dependency coefficients for each dependent reaction by reaction ID.
    """
    # SECTION: identify deterministic basis and dependent columns
    independent_indexes = independent_reaction_indices(matrix)
    dependent_indexes = dependent_reaction_indices(
        matrix=matrix,
        independent_indices_value=independent_indexes,
    )

    dependencies: ReactionDependencyMap = {}

    # SECTION: solve each dependent reaction against the same basis
    for dependent_index in dependent_indexes:
        coefficients = reaction_dependency_coefficients_by_index(
            matrix=matrix,
            independent_indices_value=independent_indexes,
            dependent_index=dependent_index,
            atol=atol,
            rtol=rtol,
        )
        dependencies[ids[dependent_index]] = {
            ids[index]: coefficient
            for index, coefficient in coefficients.items()
        }

    return dependencies


def independent_reaction_ids(
    matrix: StoichiometricMatrix,
    ids: list[str],
) -> list[str]:
    """
    Return a deterministic independent reaction basis by incremental rank test.
    """
    # NOTE: names are derived from the deterministic index basis.
    indexes = independent_reaction_indices(matrix)
    for index in indexes:
        logger.debug("Reaction %s increases network rank.", ids[index])

    return [ids[index] for index in indexes]


def dependent_reaction_ids(
    matrix: StoichiometricMatrix,
    ids: list[str],
) -> list[str]:
    """
    Return reactions that do not increase the deterministic independent basis.
    """
    # NOTE: names are derived from deterministic dependent indices.
    indexes = dependent_reaction_indices(
        matrix=matrix,
        independent_indices_value=independent_reaction_indices(matrix),
    )
    return [ids[index] for index in indexes]


def participation_matrix(matrix: StoichiometricMatrix) -> ParticipationMatrix:
    """
    Build the binary species-reaction participation matrix.
    """
    # NOTE: nonzero stoichiometric coefficients indicate participation.
    return [[1 if value != 0 else 0 for value in row] for row in matrix]
