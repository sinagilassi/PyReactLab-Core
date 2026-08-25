from __future__ import annotations

# import libs
from typing import Any


# SECTION: reaction mode aggregation
def reaction_modes(reactions: list[Any]) -> dict[str, str]:
    """
    Return normalized reaction mode/type by reaction ID.
    """
    # NOTE: Reaction owns mode normalization; the network only aggregates it.
    return {reaction.name: reaction.reaction_type for reaction in reactions}


def reaction_mode_count(reactions: list[Any]) -> dict[str, int]:
    """
    Count normalized reaction modes across the network.
    """
    # SECTION: initialize public mode buckets
    counts = {
        "equilibrium": 0,
        "reversible": 0,
        "irreversible": 0,
    }

    # SECTION: aggregate mode counts
    for mode in reaction_modes(reactions).values():
        counts.setdefault(mode, 0)
        counts[mode] += 1

    return counts
