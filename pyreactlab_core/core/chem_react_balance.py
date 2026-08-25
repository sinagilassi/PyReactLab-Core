# import libs
import logging
import re
from typing import Dict, List, Optional, Any
# ! locals
from ..configs.constants import PERIODIC_TABLE_ELEMENTS
from ..models.reactions import Reactant, Product


# SECTION: ReactChemBalance class
class ChemReactBalance:
    """
    Chemical reaction balance utilities.

    This class provides methods for:

    - parsing elemental composition from molecular formulas
    - counting elements in individual reaction components
    - calculating total reactant elemental composition
    - calculating total product elemental composition
    - checking elemental conservation
    - checking charge conservation
    - checking overall chemical balance

    Supported formula examples
    --------------------------
    H2O
    CO2
    CH3OH
    Fe(OH)3
    Ca3(PO4)2
    Al2(SO4)3
    CuSO4*5H2O
    CuSO4·5H2O
    Fe{3+}
    SO4{2-}
    OH{-}
    e{-}
    """

    # NOTE: use a set for fast membership checks
    _periodic_table_elements = frozenset(PERIODIC_TABLE_ELEMENTS)

    # NOTE: numerical tolerance for balance validation
    _balance_tolerance: float = 1e-12

    # ! ::: parse elemental composition
    def parse_elemental_composition(
        self,
        molecule: str
    ) -> Dict[str, int]:
        """
        Parse a molecular formula into its elemental composition.

        Parameters
        ----------
        molecule : str
            Molecular formula.

        Returns
        -------
        Dict[str, int]
            Mapping of element symbols to atom counts.

        Examples
        --------
        H2O
            -> {"H": 2, "O": 1}

        CO2
            -> {"C": 1, "O": 2}

        Fe(OH)3
            -> {"Fe": 1, "O": 3, "H": 3}

        Ca3(PO4)2
            -> {"Ca": 3, "P": 2, "O": 8}

        CuSO4*5H2O
            -> {"Cu": 1, "S": 1, "O": 9, "H": 10}

        SO4{2-}
            -> {"S": 1, "O": 4}

        e{-}
            -> {}
        """
        try:
            # SECTION: validate input
            if not isinstance(molecule, str):
                raise ValueError("Molecule must be a string.")

            molecule = molecule.strip()

            if not molecule:
                raise ValueError("Molecule cannot be empty.")

            # SECTION: electron
            # NOTE: electron contains no atomic elements
            if molecule == "e" or molecule.startswith("e{"):
                return {}

            # SECTION: remove ionic charge
            # Fe{3+} -> Fe
            # SO4{2-} -> SO4
            # OH{-} -> OH
            molecule = re.sub(
                r"\{(?:\d+)?[+-]\}$",
                "",
                molecule
            )

            # SECTION: hydrate / adduct sections
            # supports:
            # CuSO4*5H2O
            # CuSO4·5H2O
            sections = re.split(r"[*·]", molecule)

            total_composition: Dict[str, int] = {}

            for section in sections:
                section = section.strip()

                if not section:
                    continue

                # SECTION: hydrate multiplier
                # 5H2O -> multiplier=5, formula=H2O
                match = re.fullmatch(
                    r"(\d+)(.+)",
                    section
                )

                if match:
                    section_multiplier = int(match.group(1))
                    section_formula = match.group(2)
                else:
                    section_multiplier = 1
                    section_formula = section

                # NOTE: parse formula section
                section_composition = self._parse_formula_section(
                    section_formula
                )

                # NOTE: accumulate composition
                for element, count in section_composition.items():
                    total_composition[element] = (
                        total_composition.get(element, 0)
                        + count * section_multiplier
                    )

            return total_composition

        except Exception as e:
            raise Exception(
                f"Error parsing elemental composition "
                f"for molecule '{molecule}': {e}"
            )

    # ! ::: parse formula section
    def _parse_formula_section(
        self,
        formula: str
    ) -> Dict[str, int]:
        """
        Parse one molecular formula section.

        Supports:
        - chemical element symbols
        - numerical subscripts
        - parentheses
        - nested parentheses
        """

        # SECTION: parse number
        def parse_number(
            index: int
        ) -> tuple[int, int]:

            start = index

            while (
                index < len(formula)
                and formula[index].isdigit()
            ):
                index += 1

            # NOTE: implicit multiplier = 1
            if start == index:
                return 1, index

            return int(formula[start:index]), index

        # SECTION: recursive group parser
        def parse_group(
            index: int,
            nested: bool = False
        ) -> tuple[Dict[str, int], int]:

            composition: Dict[str, int] = {}

            while index < len(formula):

                char = formula[index]

                # SECTION: nested parentheses
                if char == "(":
                    group_composition, index = parse_group(
                        index + 1,
                        nested=True
                    )

                    # NOTE: group multiplier
                    multiplier, index = parse_number(index)

                    for element, count in group_composition.items():
                        composition[element] = (
                            composition.get(element, 0)
                            + count * multiplier
                        )

                    continue

                # SECTION: closing parentheses
                if char == ")":
                    if not nested:
                        raise ValueError(
                            f"Unexpected ')' in formula '{formula}'."
                        )

                    return composition, index + 1

                # SECTION: element
                if char.isupper():
                    element = char
                    index += 1

                    # NOTE: optional lowercase character
                    if (
                        index < len(formula)
                        and formula[index].islower()
                    ):
                        element += formula[index]
                        index += 1

                    # SECTION: validate element
                    if element not in self._periodic_table_elements:
                        raise ValueError(
                            f"Unknown chemical element '{element}' "
                            f"in formula '{formula}'."
                        )

                    # NOTE: element subscript
                    multiplier, index = parse_number(index)

                    composition[element] = (
                        composition.get(element, 0)
                        + multiplier
                    )

                    continue

                # SECTION: invalid character
                raise ValueError(
                    f"Unexpected character '{char}' "
                    f"in formula '{formula}'."
                )

            # NOTE: nested group must have closing parenthesis
            if nested:
                raise ValueError(
                    f"Unclosed '(' in formula '{formula}'."
                )

            return composition, index

        # SECTION: parse formula
        composition, final_index = parse_group(0)

        if final_index != len(formula):
            raise ValueError(
                f"Could not fully parse formula '{formula}'."
            )

        return composition

    # ! ::: count elements
    def count_elements(
        self,
        molecule: str,
        coefficient: float = 1.0
    ) -> Dict[str, float]:
        """
        Count atoms in a molecule including its stoichiometric coefficient.

        Parameters
        ----------
        molecule : str
            Molecular formula.

        coefficient : float
            Stoichiometric coefficient.

        Returns
        -------
        Dict[str, float]
            Element counts multiplied by the stoichiometric coefficient.

        Examples
        --------
        count_elements("H2O", 2)

        returns

        {
            "H": 4.0,
            "O": 2.0
        }
        """
        try:
            # SECTION: validate coefficient
            if not isinstance(coefficient, (int, float)):
                raise ValueError(
                    "Coefficient must be an integer or float."
                )

            if coefficient < 0:
                raise ValueError(
                    "Coefficient cannot be negative."
                )

            # SECTION: molecular composition
            composition = self.parse_elemental_composition(
                molecule
            )

            # NOTE: multiply by stoichiometric coefficient
            return {
                element: count * coefficient
                for element, count in composition.items()
            }

        except Exception as e:
            raise Exception(
                f"Error counting elements in molecule "
                f"'{molecule}': {e}"
            )

    # ! ::: count reaction side elements
    def count_side_elements(
        self,
        components: List[Reactant] | List[Product]
    ) -> Dict[str, float]:
        """
        Calculate total elemental composition of one reaction side.

        Parameters
        ----------
        components : List[Reactant] | List[Product]
            Reactants or products.

        Returns
        -------
        Dict[str, float]
            Total element counts.
        """
        try:
            totals: Dict[str, float] = {}

            for component in components:

                composition = self.count_elements(
                    molecule=component["molecule"],
                    coefficient=component["coefficient"]
                )

                for element, count in composition.items():
                    totals[element] = (
                        totals.get(element, 0.0)
                        + count
                    )

            return totals

        except Exception as e:
            raise Exception(
                f"Error counting reaction-side elements: {e}"
            )

    # ! ::: elemental balance
    def elemental_balance(
        self,
        reactants: List[Reactant],
        products: List[Product]
    ) -> Dict[str, Any]:
        """
        Calculate elemental conservation across a reaction.

        Returns
        -------
        Dict[str, Any]

        Example
        -------
        {
            "reactant_elements": {
                "C": 1,
                "O": 2,
                "H": 6
            },

            "product_elements": {
                "C": 1,
                "O": 2,
                "H": 6
            },

            "net_elements": {
                "C": 0,
                "H": 0,
                "O": 0
            },

            "is_element_balanced": True
        }
        """
        try:
            # SECTION: reactant elemental totals
            reactant_elements = self.count_side_elements(
                reactants
            )

            # SECTION: product elemental totals
            product_elements = self.count_side_elements(
                products
            )

            # SECTION: all participating elements
            all_elements = (
                set(reactant_elements.keys())
                | set(product_elements.keys())
            )

            # SECTION: net elemental balance
            # product - reactant
            net_elements = {
                element: (
                    product_elements.get(element, 0.0)
                    - reactant_elements.get(element, 0.0)
                )
                for element in sorted(all_elements)
            }

            # SECTION: balance check
            is_element_balanced = all(
                abs(value) <= self._balance_tolerance
                for value in net_elements.values()
            )

            return {
                "reactant_elements": reactant_elements,
                "product_elements": product_elements,
                "net_elements": net_elements,
                "is_element_balanced": is_element_balanced,
            }

        except Exception as e:
            raise Exception(
                f"Error calculating elemental balance: {e}"
            )

    # ! ::: count side charge
    def count_side_charge(
        self,
        components: List[Reactant] | List[Product]
    ) -> float:
        """
        Calculate total electrical charge for one reaction side.
        """
        try:
            total_charge = sum(
                component["charge"]
                * component["coefficient"]
                for component in components
            )

            return total_charge

        except Exception as e:
            raise Exception(
                f"Error calculating reaction-side charge: {e}"
            )

    # ! ::: charge balance
    def charge_balance(
        self,
        reactants: List[Reactant],
        products: List[Product]
    ) -> Dict[str, Any]:
        """
        Calculate electrical charge conservation.

        Returns
        -------
        Dict[str, Any]

        Example
        -------
        {
            "total_reactant_charge": 0,
            "total_product_charge": 0,
            "net_charge": 0,
            "is_charge_balanced": True
        }
        """
        try:
            # SECTION: charge totals
            total_reactant_charge = self.count_side_charge(
                reactants
            )

            total_product_charge = self.count_side_charge(
                products
            )

            # NOTE: products - reactants
            net_charge = (
                total_product_charge
                - total_reactant_charge
            )

            is_charge_balanced = (
                abs(net_charge)
                <= self._balance_tolerance
            )

            return {
                "total_reactant_charge": total_reactant_charge,
                "total_product_charge": total_product_charge,
                "net_charge": net_charge,
                "is_charge_balanced": is_charge_balanced,
            }

        except Exception as e:
            raise Exception(
                f"Error calculating charge balance: {e}"
            )

    # ! ::: reaction balance
    def reaction_balance(
        self,
        reactants: List[Reactant],
        products: List[Product]
    ) -> Dict[str, Any]:
        """
        Perform complete chemical reaction balance validation.

        Checks:
        - elemental conservation
        - electrical charge conservation

        Returns
        -------
        Dict[str, Any]
        """
        try:
            # SECTION: elemental balance
            element_result = self.elemental_balance(
                reactants=reactants,
                products=products
            )

            # SECTION: charge balance
            charge_result = self.charge_balance(
                reactants=reactants,
                products=products
            )

            # SECTION: overall balance
            is_balanced = (
                element_result["is_element_balanced"]
                and charge_result["is_charge_balanced"]
            )

            return {
                **element_result,
                **charge_result,
                "is_balanced": is_balanced,
            }

        except Exception as e:
            raise Exception(
                f"Error validating chemical reaction balance: {e}"
            )
