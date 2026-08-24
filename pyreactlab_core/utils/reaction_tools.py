# import libs
import logging
from typing import Any, Dict, List, Optional, cast, NamedTuple
# locals
from ..configs.constants import (
    IRREVERSIBLE_REACTION_MODE_SYMBOLS,
    REVERSIBLE_REACTION_MODE_SYMBOLS,
    EQUILIBRIUM_REACTION_MODE_SYMBOLS,
    ReactionMode,
    ReactionDirection,
    ReactionType,
    REACTION_SYMBOLIC_MODES,
)

# NOTE: logger
logger = logging.getLogger(__name__)

# ! reaction mode


class ReactionSymbolInfo(NamedTuple):
    reaction: str
    symbol: str
    type: ReactionType
    mode: ReactionMode
    direction: ReactionDirection


def check_reaction(
    reaction: str,
) -> ReactionSymbolInfo:
    """
Determine the reaction mode symbol, type, and direction from a given reaction string.

    Parameters
    ----------
    reaction : str
        The chemical reaction equation as a string.

    Returns
    -------
    ReactionSymbolInfo
        A named tuple containing the reaction mode symbol, type, and direction.
    """
    try:
        # SECTION: check for reaction mode symbols in the reaction string
        rxn_direction = None
        rxn_mode = None
        rxn_symbol = None

        # Longest symbols first to avoid "=" matching "<=>"
        symbols = sorted(
            REACTION_SYMBOLIC_MODES,
            key=len,
            reverse=True,
        )

        for symbol in symbols:
            if symbol in reaction:
                mode, direction = REACTION_SYMBOLIC_MODES[symbol]
                rxn_mode = mode
                rxn_direction = direction
                rxn_symbol = symbol
                # replace the symbol with :::
                reaction = reaction.replace(symbol, ":::")
                # break after finding the first valid symbol
                break

        # NOTE: if a valid reaction mode symbol is found, return it
        if (
            rxn_symbol is not None and
            rxn_symbol is not None and
            rxn_direction is not None
        ):
            if rxn_mode == "irreversible":
                return ReactionSymbolInfo(
                    reaction=reaction.replace(":::", '=>'),
                    symbol=rxn_symbol,
                    type=cast(ReactionType, "irreversible"),
                    mode=cast(ReactionMode, "=>"),
                    direction=cast(ReactionDirection, "forward")
                )
            elif rxn_mode == "reversible":
                return ReactionSymbolInfo(
                    reaction=reaction.replace(":::", '<=>'),
                    symbol=rxn_symbol,
                    type=cast(ReactionType, "reversible"),
                    mode=cast(ReactionMode, "<=>"),
                    direction=cast(ReactionDirection, "forward")
                )
            elif rxn_mode == "equilibrium":
                return ReactionSymbolInfo(
                    reaction=reaction.replace(":::", '='),
                    symbol=rxn_symbol,
                    type=cast(ReactionType, "equilibrium"),
                    mode=cast(ReactionMode, "="),
                    direction=cast(ReactionDirection, "forward")
                )
            else:
                raise ValueError(f"Unknown reaction type: {rxn_mode}")

        # NOTE: no valid reaction mode found
        raise
    except Exception as e:
        raise Exception(f"Error determining reaction mode: {e}")

# ! get reaction mode symbol


def get_reaction_mode_symbol(reaction: str) -> ReactionMode:
    """
    Determine the reaction mode symbol from a given reaction string.

    Parameters
    ----------
    reaction : str
        The chemical reaction equation as a string.

    Returns
    -------
    ReactionMode
        The reaction mode symbol.
    """
    try:
        # SECTION: check for reaction mode symbols in the reaction string
        rxn_symbol_info = check_reaction(reaction)
        return rxn_symbol_info.mode
    except Exception as e:
        raise Exception(f"Error determining reaction mode symbol: {e}")

# ! get reaction expression


def normalize_reaction_expression(reaction: str) -> str:
    """
    Normalize the reaction expression by removing the reaction mode symbol and updating the symbol to a standard form including "<=>", "=>", or "=".

    Parameters
    ----------
    reaction : str
        The chemical reaction equation as a string.

    Returns
    -------
    str
        The reaction expression without the mode symbol.
    """
    try:
        # SECTION: check for reaction mode symbols in the reaction string
        rxn_symbol_info = check_reaction(reaction)
        return rxn_symbol_info.reaction
    except Exception as e:
        raise Exception(f"Error determining reaction expression: {e}")
