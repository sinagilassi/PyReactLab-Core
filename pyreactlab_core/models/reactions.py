# import libs
from __future__ import annotations
from typing import TypedDict, Literal, Dict, List, Optional
from pydantic import BaseModel, Field
from pythermodb_settings.models import Component
# locals

# SECTION: Models
# NOTE: Phase Rule
PhaseRule = Literal["gas", "liquid", "aqueous", "solid"]

# NOTE: reactants


class Reactant(TypedDict):
    coefficient: float
    molecule: str
    charge: int
    state: str
    molecule_state: str

# NOTE: products


class Product(TypedDict):
    coefficient: float
    molecule: str
    charge: int
    state: str
    molecule_state: str
