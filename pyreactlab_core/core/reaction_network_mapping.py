from __future__ import annotations

# import libs
from typing import Any


# SECTION: component and reaction mapping
def reaction_component_map(
    reactions: list[Any],
    component_ids: list[str],
) -> dict[str, list[str]]:
    """
    Return participating component IDs by reaction ID.
    """
    # NOTE: set membership is internal only; output order follows reaction parsing.
    component_set = set(component_ids)
    mapping: dict[str, list[str]] = {}

    # SECTION: map each reaction to its participating species
    for reaction in reactions:
        mapping[reaction.name] = [
            component_id
            for component_id in reaction.reactants_names + reaction.products_names
            if component_id in component_set
        ]

    return mapping


def component_reaction_map(
    reactions: list[Any],
    component_ids: list[str],
) -> dict[str, list[str]]:
    """
    Return reaction IDs by participating component ID.
    """
    # SECTION: invert the reaction-to-component map
    reaction_map = reaction_component_map(reactions, component_ids)
    mapping = {component_id: [] for component_id in component_ids}

    for reaction_id, reaction_components in reaction_map.items():
        for component_id in reaction_components:
            mapping[component_id].append(reaction_id)

    return mapping


def component_metadata(reactions: list[Any]) -> dict[str, dict[str, Any]]:
    """
    Return first-seen component metadata keyed by molecule-state ID.
    """
    # NOTE: metadata is sourced from Reaction reactant/product dictionaries.
    metadata: dict[str, dict[str, Any]] = {}

    # SECTION: collect unique component metadata
    for reaction in reactions:
        for item in reaction.reactants + reaction.products:
            component_id = item["molecule_state"]
            if component_id not in metadata:
                metadata[component_id] = dict(item)

    return metadata
