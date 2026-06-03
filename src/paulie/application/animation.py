"""
    Module for animating the transformation of the anti-commutation graph into a canonical form.
"""
from matplotlib.animation import Animation
from paulie.helpers._recording import RecordGraph
from paulie.helpers.drawing import _animation_graph
from paulie.common.pauli_string_collection import PauliStringCollection


def animation_anti_commutation_graph(
    generators: PauliStringCollection,
    storage: dict[str, str] | None = None,
    interval: int = 1000,
    show: bool = False,
) -> Animation:
    """
    Generates an animation showing the transformation of the anti-commutation
    graph into canonical form.

    The animation is driven by a
    :class:`~paulie.classifier.recording_canonicalizer.RecordingCanonicalizer`, which observes the
    classification algorithm and records each step. Use the colour legend in
    :data:`paulie.helpers.drawing.NODE_ROLE_COLORS` to interpret node roles.

    Args:
        generators (PauliStringCollection): Collection of Pauli strings.
        storage (dict[str, str], optional): Location and format to save the
            animation to. Expected keys are ``"filename"`` and ``"writer"``.
        interval (int, optional): Interval between recording frames in
            milliseconds.
        show (bool, optional): Whether to display the animation window.

    Returns:
        matplotlib.animation.Animation
    """
    record = RecordGraph()
    generators.set_record(record)
    generators.classify()
    generators.set_record(None)
    return _animation_graph(
        record,
        interval=interval,
        storage=storage,
        show=show,
    )
