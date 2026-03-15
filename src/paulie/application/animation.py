"""
    Module for animating the transformation of the anti-commutation graph into a canonical form.
"""
from matplotlib.animation import Animation
from paulie.helpers._recording import RecordGraph
from paulie.helpers.drawing import _animation_graph
from paulie.common.pauli_string_collection import PauliStringCollection



def animation_anti_commutation_graph(
    generators: PauliStringCollection,
    storage=None,
    interval: int = 1000,
    show: bool = False,
) -> Animation:
    """
    Generates an animation showing the transformation of the anti-commutation
    graph into canonical form.

    Args:
        generators (PauliStringCollection): Collection of Pauli strings.
        storage (dict[str, str], optional): Location and format to save the
            animation to.
        interval (int, optional): Interval between recording frames in
            milliseconds.
        show (bool, optional): Whether to display the animation window.

    Returns:
        matplotlib.animation.Animation
    """
    record = RecordGraph()
    generators.set_record(record)
    generators.get_class()
    return _animation_graph(
        record,
        interval=interval,
        storage=storage,
        show=show,
    )
