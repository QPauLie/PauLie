.. module:: paulie.classifier

Classifier (:mod:`paulie.classifier`)
=====================================

Algebra classification
----------------------
Utilities for classifying the DLA of a set of Pauli strings.

.. autosummary::
   :toctree: generated/

   classification.Classification
   classification.Morph
   classification.TypeAlgebra
   classification.TypeGraph
   morph_factory.MorphFactory
   recording_morph_factory.RecordingMorphFactory

Exceptions
----------
Exceptions thrown by the above methods.

.. autosummary::
   :toctree: generated/

   classification.ClassificationException
   morph_factory.MorphFactoryException
   morph_factory.NotConnectedException
   morph_factory.RaiseException
   morph_factory.AppendedException
   morph_factory.CheckAppendedException
   morph_factory.DependentException