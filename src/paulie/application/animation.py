"""
Animation building transformation anti-commutation graph
"""
from paulie.helpers.recording import RecordGraph
from paulie.helpers.drawing import animation_graph
from paulie.common.pauli_string_collection import PauliStringCollection


def animation_anti_commutation_graph(generators: PauliStringCollection,
                                     storage=None, interval:int =1000) -> None:
    """
    Generates an animation showing the transformation of the anti-commutation graph
    into canonical form.

    Args:
        generators (PauliStringCollection): Collection of Pauli strings.
        storage (dict): Storage in file.

            - filename (string): path to file
            - writer (string): Specifies the software used to write the animation.
              Common options include `ffmpeg` for MP4, AVI, etc. output, and `pillow`
              and `imagemagick` for GIF output. Note that the necessary libraries must
              be available and usable by `matplotlib`.

        interval (int): Interval between recording frames.

    Returns:
        None
    """
    record = RecordGraph()
    generators.set_record(record)
    generators.get_class()
    animation_graph(record, interval=interval, storage=storage)
