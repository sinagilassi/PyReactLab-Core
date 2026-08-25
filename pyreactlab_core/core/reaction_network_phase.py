from __future__ import annotations

# import libs
from typing import Any

# locals
from .reaction_network_mapping import component_metadata
from .reaction_network_species import unique_preserve_order


# SECTION: phase analysis
def phases(reactions: list[Any]) -> list[str]:
    """
    Return unique phase/state labels in first-seen component order.
    """
    # SECTION: collect phases from component metadata
    values: list[str] = []
    for data in component_metadata(reactions).values():
        values.append(data["state"])

    return unique_preserve_order(values)


def components_by_phase(
    reactions: list[Any],
    component_ids: list[str],
) -> dict[str, list[str]]:
    """
    Return component IDs grouped by phase/state.
    """
    # NOTE: component order inside each phase follows component_ids.
    metadata = component_metadata(reactions)
    mapping = {phase: [] for phase in phases(reactions)}

    # SECTION: group components by parsed state
    for component_id in component_ids:
        phase = metadata[component_id]["state"]
        mapping.setdefault(phase, []).append(component_id)

    return mapping


def reactions_by_phase(reactions: list[Any]) -> dict[str, list[str]]:
    """
    Return reaction IDs grouped by phases that participate in each reaction.
    """
    # NOTE: reactions can appear in more than one phase bucket.
    mapping = {phase: [] for phase in phases(reactions)}

    # SECTION: collect phase memberships per reaction
    for reaction in reactions:
        reaction_phases = unique_preserve_order(
            [
                item["state"]
                for item in reaction.reactants + reaction.products
            ]
        )
        for phase in reaction_phases:
            mapping.setdefault(phase, []).append(reaction.name)

    return mapping


def single_phase_reactions(reactions: list[Any]) -> list[str]:
    """
    Return reaction IDs whose species all share one phase/state.
    """
    # SECTION: classify reactions by unique participating phases
    values: list[str] = []
    for reaction in reactions:
        reaction_phases = {
            item["state"]
            for item in reaction.reactants + reaction.products
        }
        if len(reaction_phases) == 1:
            values.append(reaction.name)

    return values


def multiphase_reactions(reactions: list[Any]) -> list[str]:
    """
    Return reaction IDs whose species span more than one phase/state.
    """
    # SECTION: classify reactions by unique participating phases
    values: list[str] = []
    for reaction in reactions:
        reaction_phases = {
            item["state"]
            for item in reaction.reactants + reaction.products
        }
        if len(reaction_phases) > 1:
            values.append(reaction.name)

    return values
