"""
Animation building transformation anti-commutation graph
"""
from paulie.helpers.recording import RecordGraph
from paulie.helpers.drawing import animation_graph
from paulie.common.pauli_string_collection import PauliStringCollection


def animation_anti_commutation_graph(generators: PauliStringCollection,
                                     storage=None, interval:int =1000) -> None:
    """
    Animation building transformation anti-commutation graph
    generators - list of generators
    Args:
        generators: collection of Pauli strings
        storage: storage in file: 
                 dictionary: filename: path to file
                             writer: specifies the software used
                                     to write the animation. 
                                     Common options include:
                                    'ffmpeg' (for MP4, AVI, etc.,
                                     requires FFmpeg to be installed and accessible
                                     in your system's PATH, or its path specified via
                                     matplotlib.rcParams['animation.ffmpeg_path']).
                                    'pillow' (for GIF, requires the Pillow library).
                                    'imagemagick' (for GIF, requires ImageMagick).

        interval: interval between recording frames
    Returns:
        None
    """
    record = RecordGraph()
    generators.set_record(record)
    generators.get_class()
    animation_graph(record, interval=interval, storage=storage)
