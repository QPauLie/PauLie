"""
    Module for animating the transformation of the anti-commutation graph into a canonical form.
"""
from paulie.helpers._recording import RecordGraph
from paulie.helpers.drawing import _animation_graph
from paulie.common.pauli_string_collection import PauliStringCollection


def animation_anti_commutation_graph(generators: PauliStringCollection,
                                     storage=None, interval:int =1000) -> None:
    """
    Generates an animation showing the transformation of the anti-commutation graph into canonical
    form.

    Args:
        generators (PauliStringCollection): Collection of Pauli strings.
        storage (dict, optional): Location and format to save the animation to. Defaults to `None`,
            in which case the animation is not saved.

            - filename (string): Path to the output file.
            - writer (string): Specifies the software used to write the animation. Common options
              include `ffmpeg` for MP4, AVI, etc. output, and `pillow` and `imagemagick` for GIF
              output. Note that the necessary libraries must be available and usable by
              `matplotlib`.
        interval (int, optional): Interval between recording frames in milliseconds. Defaults to
            1000 milliseconds.

    Returns:
        None
    """
    record = RecordGraph()
    generators.set_record(record)
    generators.get_class()
    _animation_graph(record, interval=interval, storage=storage)
