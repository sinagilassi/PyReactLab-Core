# import libs
from typing import List
from pythermodb_settings.models import Component, ComponentKey
from pythermodb_settings.utils import set_component_id


def map_component_ids(
        component_ids: list[str],
        component: List[Component],
        component_key: ComponentKey
):
    """
    Map component IDs to their corresponding components.

    Parameters
    ----------
    component_ids : list[str]
        A list of component IDs to be mapped.
    component : List[Component]
        A list of Component objects to be mapped against.
    component_key : ComponentKey
        The key used to identify components in the reaction.

    Returns
    -------
    Dict[str, Component]
        A dictionary mapping component IDs to their corresponding Component objects.
    """
    id_to_component = {}

    # iterate over components and map ids
    for comp in component:
        comp_id = set_component_id(comp, component_key)
        if comp_id in component_ids:
            id_to_component[comp_id] = comp

    return id_to_component
