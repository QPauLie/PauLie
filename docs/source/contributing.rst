Contribute to PauLie
####################

We welcome contributions from everyone! By participating in PauLie, you agree to follow our guidelines.

1. As stated in `Getting started <getting_started.html>`_, ensure you have Python 3.12 or higher installed on your machine or in a virtual environment.
2. You requires to `install uv <https://docs.astral.sh/uv/getting-started/installation/#pypi>`_ to manage the dependencies of PauLie.
3. Make sure you have a `GitHub account <https://github.com/signup/free>`_.
4. `Fork <https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/working-with-forks/fork-a-repo>`_ the repository on GitHub.
5. `Clone <https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository>`_ your fork locally:

.. code-block:: bash

       git clone https://github.com/QPauLie/PauLie.git

6. Navigate to your local clone PauLie reprository:

.. code-block:: bash

       cd PauLie

7. Install dependencies with uv:

.. code-block:: bash

       uv sync --all-extras --dev

uv will create a virtual environment and install all dependencies there. You are now ready to run, test, and develop the project.

Making changes
==============

1. Before starting to create your ideas, make sure that your fork in up to date with the original PauLie repository:

.. code-block:: bash

    git fetch upstream
    git checkout main
    git merge upstream/main

2. Create your feature or idea on your local fork. It’s best to make changes on a dedicated `branch <https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-and-deleting-branches-within-your-repository>`_, with the branch name reflecting the feature or fix you are adding.

3. Include tests and documentation for new features. If you’re adding a new feature, ensure you provide test cases and update the documentation. See :ref:`Adding a New Feature <adding-new-features>` for a detailed checklist.

4. When the code is ready to go, make sure you run the :ref:`test suite <functional-testing>` using :code:`pytest`, :code:`ruff`, etc.

5. Open a `pull request (PR) <https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/about-pull-requests>`_ on your local fork. When your feature is ready for review, go to your fork of PauLie on GitHub and open a PR. Any subsequent commits to that branch will automatically be added to the open PR. This will commit your feature to the branch of your local fork.

6. Open a pull request on original fork. Once your changes are ready for merging, go to the `original PauLie <https://github.com/QPauLie/PauLie/pulls>`_ on GitHub. Then, open PR, selecting "compare across fork" to submit your changes from your fork’s branch, and add a comment indicating that your code is ready for review. The code will be reviewed after the continuous integration (CI) workflow passes and the primary developer approves the review.

.. _adding-new-features:

Adding New Features
===================
The new feature should follow these guidelines:

1. The function docstring should follow the Google style as outlined in `References in Docstrings <https://www.sphinx-doc.org/en/master/usage/extensions/example_google.html>`_.

2. The docstring for a new feature must include:

   - A theoretical description of the feature.
   - One or more examples in an **Examples** subsection.
   - A **References** subsection.

3. All added lines must be covered by tests and appear in the :code:`pytest` code coverage report. See :ref:`Testing <functional-testing>` for more details.

4. The code and unit tests for the new feature should follow the Code Style guidelines, using :code:`ruff` for :ref:`lint checking <lint-checking>`.

5. The new feature must be added to the __init__.py file of its module to avoid import issues.

6. If the new feature introduces a new module, it must be listed in docs/autoapi_members.rst so that it appears on the API Reference page via sphinx-autoapi.

Before submitting a pull request, please make sure your changes pass all tests.

.. _functional-testing:

Functional Testing
==================

Ensure that both :code:`pytest` and :code:`pytest-cov` are installed on your machine.
If they are not installed, follow the installation instructions for `pytest <https://docs.pytest.org/en/stable/getting-started.html>`_ and `pytest-cov <https://pytest-cov.readthedocs.io/en/latest/readme.html#installation>`_.

To verify that your code is working correctly, run the following commands in the PauLie directory:

.. code-block:: bash

    uv run pytest -q

The pytest module is used for running tests, and pytest-cov can generate local code coverage reports.
Code coverage can be checked using the :code:`--cov-report=term-missing` option in :code:`pytest-cov`.

Any changes to an existing module should have corresponding tests placed in its :code:`tests` directory (e.g., :code:`PauLie/tests/test_module.py`).

.. _lint-checking:

Lint Checking
=============

Formatting and linting are handled using :code:`ruff`. To install :code:`ruff`, refer to the official `ruff installation instructions <https://docs.astral.sh/ruff/installation/>`_.

To detect errors and style issues such as unused imports, undefined variables, style violations, and complexity problems, run:

.. code-block:: bash

    uv run ruff check

To verify that your code is formatted according to Ruff’s formatter rules, run:

.. code-block:: bash

    uv run ruff format --check

Note: without the :code:`--check` option, :code:`ruff format` will modify your code automatically.

For deeper static analysis and design-level checks, we use :code:`pylint`.

To lint the source code, test scripts, and examples without installing :code:`pylint`, run:

.. code-block:: bash

    uvx pylint "src/**/*.py" "tests/**/*.py" "examples/**/*.py"

.. _doc-building:

Document Building
=================

The sphinx is used in PauLie to build the HTML documentation, use:

.. code-block:: bash

    uv run sphinx-build -b html -W docs/source docs/_build/html

