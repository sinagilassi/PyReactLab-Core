# import libs
import logging
from typing import Any, Dict, List, Sequence

# import packages
from pythermodb_settings.models import Component, ComponentKey
from pythermodb_settings.utils import set_component_id


# SECTION: ReactionComponentMapper class
class ReactionComponentMapper:
    """Map parsed reaction participants back to configured components."""

    # NOTE: attributes are supplied by ChemReact at runtime
    components: List[Component] | None
    component_ids: List[str]
    component_keys: List[ComponentKey]
    _component_checker: bool
    _stoichiometry_source: dict[str, Any]

    def collect_components(
        self,
        reactants: Sequence[Any],
        products: Sequence[Any],
    ) -> List[Component]:
        """
        Collect components from reactants and products.
        """
        try:
            # SECTION: initialize outputs
            components: List[Component] = []
            components_ids_: List[str] = []

            # NOTE: no configured component source means there is nothing to map
            if self.components is None:
                return components

            # SECTION: combine parsed participants
            # NOTE: Sequence[Any] keeps TypedDict inputs accepted by static checkers
            reaction_participants = list(reactants) + list(products)
            reaction_participants_num = len(reaction_participants)

            # SECTION: collect matching components
            for item in reaction_participants:
                # NOTE: molecule_state is the Formula-State id from analysis
                component_id_ = item['molecule_state']

                # ! only collect configured components and avoid duplicates
                if (
                    component_id_ in self.component_ids and
                    component_id_ not in components_ids_
                ):
                    # NOTE: keep local ids in the same order as self.components
                    components_ids_.append(component_id_)
                    index = self.component_ids.index(component_id_)
                    component = self.components[index]
                    components.append(component)

            # SECTION: component availability check
            # ? duplicate reactants/products currently count as separate participants
            if (
                len(components) != reaction_participants_num
            ):
                self._component_checker = False
            else:
                self._component_checker = True

            # NOTE: return configured Component objects found in this reaction
            return components
        except Exception as e:
            raise Exception(f"Error collecting components: {e}")

    def map_components(
            self,
            reactants: Sequence[Any],
            products: Sequence[Any],
    ) -> Dict[str, Component]:
        """
        Map components from reactants and products.
        """
        try:
            # SECTION: initialize output
            component_map: Dict[str, Component] = {}

            # NOTE: no configured component source means there is nothing to map
            if self.components is None:
                return component_map

            # SECTION: combine parsed participants
            reaction_participants = list(reactants) + list(products)

            # SECTION: map molecule-state ids to Component objects
            for item in reaction_participants:
                # NOTE: molecule_state is the Formula-State id from analysis
                component_id_ = item['molecule_state']

                # ! keep the first configured component for each molecule-state id
                if (
                    component_id_ in self.component_ids and
                    component_id_ not in component_map
                ):
                    # NOTE: component_ids mirrors self.components indexing
                    index = self.component_ids.index(component_id_)
                    component = self.components[index]
                    component_map[component_id_] = component

            # NOTE: return lookup map used by Reaction.map_components
            return component_map
        except Exception as e:
            raise Exception(f"Error mapping components: {e}")

    def reaction_stoichiometry_dict(
        self,
        reaction_stoichiometry: Dict[str, float],
        component_key: ComponentKey = "Name-Formula"
    ) -> Dict[str, float]:
        """
        Build a stoichiometry dictionary keyed by the requested component key.
        """
        # SECTION: validate component source
        if self.components is None:
            return {}

        # NOTE: reaction stoichiometry is keyed by Formula-State ids
        component_ids = list(reaction_stoichiometry.keys())

        # SECTION: map stoichiometric coefficients to requested component ids
        id_to_component = {}
        for comp in self.components:
            # NOTE: Formula-State is the source key from analyze_reaction
            comp_id = set_component_id(comp, 'Formula-State')  # type: ignore
            if comp_id in component_ids:
                # NOTE: convert to requested component key format
                comp_id_new = set_component_id(
                    comp,
                    component_key
                )
                id_to_component[comp_id_new] = reaction_stoichiometry[comp_id]
            else:
                # ! configured component is not present in this reaction
                logging.warning(
                    f"Component '{comp_id}' not found in reaction stoichiometry.")

        # NOTE: return stoichiometry keyed by requested component identifier
        return id_to_component

    def build_stoichiometry_source(
            self,
            reaction_stoichiometry: Dict[str, float],
    ) -> Dict[str, Any]:
        """
        Build the stoichiometry source for the reaction component keys.
        """
        # SECTION: validate component source
        if self.components is None:
            return {}

        # SECTION: build stoichiometry source
        stoichiometry_source = {}

        # NOTE: generate one stoichiometry dict per requested component key
        for key in self.component_keys:
            id_to_component = self.reaction_stoichiometry_dict(
                reaction_stoichiometry,
                component_key=key
            )

            stoichiometry_source[key] = id_to_component

        # NOTE: cache latest source on the owning ChemReact instance
        self._stoichiometry_source = stoichiometry_source
        return stoichiometry_source
