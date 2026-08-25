from __future__ import annotations

# import libs
from typing import Any


# SECTION: deterministic collection helpers
def unique_preserve_order(values: list[str]) -> list[str]:
    """
    Return unique values in first-seen order.
    """
    # NOTE: dict keys preserve insertion order in supported Python versions.
    return list(dict.fromkeys(values))


def reaction_ids(reactions: list[Any]) -> list[str]:
    """
    Return network reaction identifiers in user-provided order.
    """
    # NOTE: reaction names are the public reaction identifiers.
    return [reaction.name for reaction in reactions]


def reactants(reactions: list[Any]) -> list[str]:
    """
    Return unique network reactant species identifiers.
    """
    # SECTION: collect reactants from already-analyzed Reaction objects
    values: list[str] = []
    for reaction in reactions:
        values.extend(reaction.reactants_names)

    # NOTE: preserve first-seen network order.
    return unique_preserve_order(values)


def products(reactions: list[Any]) -> list[str]:
    """
    Return unique network product species identifiers.
    """
    # SECTION: collect products from already-analyzed Reaction objects
    values: list[str] = []
    for reaction in reactions:
        values.extend(reaction.products_names)

    # NOTE: preserve first-seen network order.
    return unique_preserve_order(values)


def all_species(reactions: list[Any]) -> list[str]:
    """
    Return unique network species identifiers.
    """
    # NOTE: reactant-first ordering is deterministic and user-facing.
    return unique_preserve_order(reactants(reactions) + products(reactions))


# SECTION: species role classification
def source_species(reactions: list[Any]) -> list[str]:
    """
    Return species that appear as reactants and never as products.
    """
    # NOTE: set membership is used only for filtering, not output ordering.
    product_set = set(products(reactions))
    return [species for species in reactants(reactions) if species not in product_set]


def sink_species(reactions: list[Any]) -> list[str]:
    """
    Return species that appear as products and never as reactants.
    """
    # NOTE: set membership is used only for filtering, not output ordering.
    reactant_set = set(reactants(reactions))
    return [species for species in products(reactions) if species not in reactant_set]


def intermediate_species(reactions: list[Any]) -> list[str]:
    """
    Return species that appear as both reactants and products.
    """
    # NOTE: keep reactant-side order for deterministic intermediate ordering.
    product_set = set(products(reactions))
    return [species for species in reactants(reactions) if species in product_set]
