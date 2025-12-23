# import libs
import logging
from typing import Literal, Optional

# setup logger
logger = logging.getLogger(__name__)

# NOTE: check if value is number


def is_number(value: str) -> bool:
    try:
        float(value)
        return True
    except (ValueError, TypeError):
        return False

# NOTE: check if value is integer


def is_integer(value: str) -> bool:
    try:
        num = float(value)   # handles "10", "10.5", "1e-3", etc.
    except (ValueError, TypeError):
        return False

    if num.is_integer():      # built-in float method
        return True
    else:
        return False
