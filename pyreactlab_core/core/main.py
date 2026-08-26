# import libs
import logging
from typing import Dict, Optional
# ! locals
from .chem_react_balance import ChemReactBalance

# NOTE: logger
logger = logging.getLogger(__name__)

# ! ::: Parse Elemental Composition


def parse_elemental_composition(
        compound: str
) -> Optional[Dict[str, int]]:
    try:
        # SECTION: validate input
        if not isinstance(compound, str) or not compound:
            logger.error(
                f"Invalid compound input: '{compound}'. Must be a non-empty string.")
            return None

        # SECTION: init
        chem_react_balance = ChemReactBalance()

        # SECTION: parse elemental composition
        element_composition = chem_react_balance.parse_elemental_composition(
            compound
        )

        return element_composition
    except Exception as e:
        logger.error(
            f"Error parsing elemental composition for compound '{compound}': {e}")
        return None


# ! ::: Parse Ionic Charge
def parse_ionic_charge(
        compound: str
) -> Optional[int]:
    try:
        # SECTION: validate input
        if not isinstance(compound, str) or not compound:
            logger.error(
                f"Invalid compound input: '{compound}'. Must be a non-empty string.")
            return None

        # SECTION: init
        chem_react_balance = ChemReactBalance()

        # SECTION: parse ionic charge
        ionic_charge = chem_react_balance.parse_ionic_charge(
            compound
        )

        return ionic_charge
    except Exception as e:
        logger.error(
            f"Error parsing ionic charge for compound '{compound}': {e}")
        return None
