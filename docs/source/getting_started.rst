.. _getting_started_reference-label:

Getting started
===============

Installing
----------

Make sure to have `python >= 3.12` installed.

Once you have `uv` installed, run:

.. code-block:: bash

   uv sync --all-extras --dev
   uv shell

Development
-----------

Typing
~~~~~~

Check types and catch errors with mypy:

.. code-block:: bash

   uv run mypy src/paulie/common src/paulie/classifier/classification.py