# import libs
import logging
from typing import List, Optional, Dict
# locals
from .models.reaction import Reaction
from .docs.chem_utils import build_stoichiometry_matrix

# NOTE: configure logger
logger = logging.getLogger(__name__)


def rxn(
        reaction_str: str,
        name: Optional[str] = None
) -> Optional[Reaction]:
    """
    Create and analyze a chemical reaction.

    Parameters
    ----------
    reaction_str : str
        reaction expression string such as "A + B => C + D"
    name : Optional[str]
        Name of the reaction

    Returns
    -------
    Reaction
        An instance of the Reaction class containing analysis results.
    """
    try:
        # SECTION: validate inputs
        if not isinstance(reaction_str, str) or not reaction_str.strip():
            logger.error("Invalid reaction string provided.")
            return None

        # check name format
        if name is not None and not isinstance(name, str):
            logger.error("Invalid name format provided.")
            return None

        # SECTION: create Reaction instance
        reaction = Reaction(
            name=name if name else "Unnamed Reaction",
            reaction=reaction_str
        )
        return reaction
    except Exception as e:
        logger.error(f"Error creating Reaction instance: {e}")
        return None


def rxn_stoichiometry(
        reaction: Reaction,
) -> Optional[List[float]]:
    """
    Get the reaction stoichiometry matrix for a given chemical reaction.

    Parameters
    ----------
    reaction : Reaction
        An instance of the Reaction class.

    Returns
    -------
    Optional[List[float]]
        The reaction stoichiometry matrix if successful, None otherwise.
    """
    try:
        # SECTION: validate inputs
        if not isinstance(reaction, Reaction):
            logger.error("Invalid Reaction instance provided.")
            return None

        # SECTION: retrieve stoichiometry matrix
        if reaction:
            return reaction.reaction_stoichiometry_matrix
        else:
            logger.error("Failed to create Reaction instance.")
            return None
    except Exception as e:
        logger.error(f"Error retrieving reaction stoichiometry matrix: {e}")
        return None


def rxns_stoichiometry(
        reactions: List[Reaction],
) -> Optional[Dict]:
    """
    Get the reaction stoichiometry matrices for a list of chemical reactions.

    Parameters
    ----------
    reactions : List[Reaction]
        List of Reaction instances.

    Returns
    -------
    Optional[List[List[float]]]
        List of reaction stoichiometry matrices if successful, None otherwise.
    """
    try:
        # SECTION: validate inputs
        if (
            not isinstance(reactions, list) or
            not all(isinstance(rxn, Reaction) for rxn in reactions)
        ):
            logger.error("Invalid reactions list provided.")
            return None

        # SECTION: retrieve stoichiometry matrices
        result = build_stoichiometry_matrix(reactions=reactions)

        # NOTE: res
        return {
            "components": result["component_list"],
            "component_ids": result["component_dict"],
            "stoichiometry_matrices_list": result["comp_coeff"],
            "stoichiometry_matrices_dict": result["comp_list"],
        }
    except Exception as e:
        logger.error(f"Error retrieving reaction stoichiometry matrices: {e}")
        return None
