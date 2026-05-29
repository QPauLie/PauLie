.. _getting_started_reference-label:

Getting started
===============

Installing
----------

`PauLie` package requires Python 3.12 or newer.
To install `PauLie`, run:

.. code-block:: bash

   pip install paulie

Development
-----------

Typing
~~~~~~

Once you have `uv` installed, set up the development environment by running:

.. code-block:: bash

   uv sync --all-extras --dev
   uv shell

Check types and catch errors with mypy:

.. code-block:: bash

   uv run mypy src/paulie/common src/paulie/classifier/classification.py