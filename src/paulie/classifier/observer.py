"""
Observer machinery for recording the canonicalization process.
"""
import enum


class CanonicalizerEvent(enum.Enum):
    """
    Event emitted by the canonicalizer at each step of the transformation.

    The value of each member is the default human-readable title of the event.
    """
    ANTICOMMUTATION_GRAPH = "Anticommutation graph"
    CENTRAL_VERTEX = "Central vertex"
    LIGHTING = "Adding"
    BUILD_CORE = "Build core"
    CONVERT_SINGLE_LIT_STATE = "Single legs in different lit states"
    TRANSFER_LIGHTNING = "Transfer lightning to long leg"
    REDUCE_LIGHTNING = "Reduce lightning on long leg"
    REMOVE_DEPENDENT = "Remove dependent vertices"
    CANONICAL = "Canonical graph"


class CanonicalizerObserver:
    """
    Base class for observers of the canonicalization process.
    """
    def update(self, event: CanonicalizerEvent, data: dict) -> None:  # pylint: disable=unused-argument
        """
        React to an event emitted by the canonicalizer.

        Args:
            event (CanonicalizerEvent): Type of event emitted.
            data (dict): Contextual data associated with the event.
        Returns:
            None
        """


class EventManager:
    """
    Subscription manager relating a canonicalizer to its observers.
    """
    def __init__(self) -> None:
        """
        Initialize an EventManager with no subscribers.

        Returns:
            None
        """
        self.subscribers: list[CanonicalizerObserver] = []

    def subscribe(self, observer: CanonicalizerObserver) -> None:
        """
        Subscribe an observer.

        Args:
            observer (CanonicalizerObserver): Observer to subscribe.
        Returns:
            None
        """
        self.subscribers.append(observer)

    def unsubscribe(self, observer: CanonicalizerObserver) -> None:
        """
        Unsubscribe an observer.

        Args:
            observer (CanonicalizerObserver): Observer to unsubscribe.
        Returns:
            None
        """
        self.subscribers.remove(observer)

    def has_subscribers(self) -> bool:
        """
        Check if there are any subscribers.

        Returns:
            bool: True if at least one observer is subscribed.
        """
        return len(self.subscribers) > 0

    def notify(self, event: CanonicalizerEvent, data: dict) -> None:
        """
        Notify all subscribers of an event.

        Args:
            event (CanonicalizerEvent): Type of event emitted.
            data (dict): Contextual data associated with the event.
        Returns:
            None
        """
        for observer in self.subscribers:
            observer.update(event, data)
