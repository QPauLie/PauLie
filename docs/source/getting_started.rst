.. _getting_started_reference-label:

===============
Getting started
===============
----------
Installing
----------
Make sure to have `python >= 3.12` installed.

Once you have `poetry` installed, run:

.. code-block:: bash

    poetry install
    poetry shell

----------
Development
----------
### Typing

Check types and catch errors with mypy:

.. code-block:: bash

    poetry run mypy src/paulie/common src/paulie/classifier/classification.py