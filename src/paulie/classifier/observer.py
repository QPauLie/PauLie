"""
    Observer pattern for the canonicalizer.

    The :class:`~paulie.classifier.canonicalizer.Canonicalizer` plays the role of the
    *publisher*. It owns an :class:`EventManager` and emits events at every relevant step of
    the classification algorithm. Objects that want to track these events (for example, the
    frame recorder used to animate the transformation) implement :class:`CanonicalizerObserver`
    and subscribe to the event manager. Subscribers only observe; they never change the
    behaviour of the algorithm.
"""
from __future__ import annotations

from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from paulie.classifier.canonicalizer import Canonicalizer


class CanonicalizerObserver:
    """
    Subscriber interface for canonicalizer events.

    Concrete subscribers override :meth:`update` to react to events emitted by a
    :class:`~paulie.classifier.canonicalizer.Canonicalizer`.
    """

    def update(self, event: str, source: "Canonicalizer", data: dict[str, Any]) -> None:
        """
        React to an event emitted by the publisher.

        Args:
            event (str): Name of the event (used as the frame title).
            source (Canonicalizer): The canonicalizer that emitted the event. Subscribers may
                read its current state (``central_vertex``, ``legs``) directly.
            data (dict[str, Any]): Contextual data describing the event, with keys matching the
                node roles tracked by the recorder (``lighting``, ``appending``, ``contracting``,
                ``p``, ``q``, ``removing``, ``replacing``, ``dependent``, ``collection``, ``init``).
        Returns:
            None
        """


class EventManager:
    """
    Subscription manager shared by composition.

    Holds the list of subscribers and forwards events to them. The canonicalizer keeps an
    instance of this class instead of inheriting the subscription logic, so the recording
    machinery can be patched in without touching the class hierarchy.
    """

    def __init__(self) -> None:
        """
        Initialize an empty subscriber list.

        Returns:
            None
        """
        self._subscribers: list[CanonicalizerObserver] = []

    def subscribe(self, observer: CanonicalizerObserver) -> None:
        """
        Register a subscriber.

        Args:
            observer (CanonicalizerObserver): Subscriber to register.
        Returns:
            None
        """
        if observer not in self._subscribers:
            self._subscribers.append(observer)

    def unsubscribe(self, observer: CanonicalizerObserver) -> None:
        """
        Remove a subscriber.

        Args:
            observer (CanonicalizerObserver): Subscriber to remove.
        Returns:
            None
        """
        if observer in self._subscribers:
            self._subscribers.remove(observer)

    def has_subscribers(self) -> bool:
        """
        Check whether any subscriber is registered.

        Returns:
            bool: True if at least one subscriber is registered.
        """
        return bool(self._subscribers)

    def notify(self, event: str, source: "Canonicalizer", data: dict[str, Any]) -> None:
        """
        Notify all subscribers of an event.

        Args:
            event (str): Name of the event.
            source (Canonicalizer): The canonicalizer emitting the event.
            data (dict[str, Any]): Contextual data describing the event.
        Returns:
            None
        """
        for observer in self._subscribers:
            observer.update(event, source, data)
