# import libs
import logging
from typing import Any, Dict, List, Optional, cast
# locals
from ..configs.constants import (
    IRREVERSIBLE_REACTION_MODE_SYMBOLS,
    REVERSIBLE_REACTION_MODE_SYMBOLS,
    EQUILIBRIUM_REACTION_MODE_SYMBOLS,
    ReactionMode,
    REACTION_SYMBOLIC_MODES,
)

# NOTE: logger
logger = logging.getLogger(__name__)

# ! reaction mode


def get_reaction_mode_symbol(
    reaction: str,
) -> ReactionMode:
    """
    Determine the reaction mode based on the reaction string.

    Parameters
    ----------
    reaction : str
        The chemical reaction equation as a string.

    Returns
    -------
    Optional[ReactionMode]
        The reaction mode symbol if found, otherwise None.
    """
    try:
        # SECTION: check for reaction mode symbols in the reaction string
        rxn_type = None

        # iterate through reversible reaction mode symbols
        for symbol in REVERSIBLE_REACTION_MODE_SYMBOLS:
            if symbol in reaction:
                # ! return the reaction mode symbol
                rxn_type = REACTION_SYMBOLIC_MODES[symbol][0]
                break  # exit loop once a match is found

        # iterate through irreversible reaction mode symbols
        if rxn_type is None:
            for symbol in IRREVERSIBLE_REACTION_MODE_SYMBOLS:
                if symbol in reaction:
                    # ! return the reaction mode symbol
                    rxn_type = REACTION_SYMBOLIC_MODES[symbol][0]
                    break  # exit loop once a match is found

        # iterate through equilibrium reaction mode symbols
        if rxn_type is None:
            for symbol in EQUILIBRIUM_REACTION_MODE_SYMBOLS:
                if symbol in reaction:
                    # ! return the reaction mode symbol
                    rxn_type = REACTION_SYMBOLIC_MODES[symbol][0]
                    break  # exit loop once a match is found

        # NOTE: if a valid reaction mode symbol is found, return it
        if rxn_type is not None:
            if rxn_type == "irreversible":
                return cast(ReactionMode, "=>")
            elif rxn_type == "reversible":
                return cast(ReactionMode, "<=>")
            elif rxn_type == "equilibrium":
                return cast(ReactionMode, "=")
            else:
                raise ValueError(f"Unknown reaction type: {rxn_type}")

        # NOTE: no valid reaction mode found
        raise
    except Exception as e:
        raise Exception(f"Error determining reaction mode: {e}")
