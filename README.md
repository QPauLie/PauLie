# PauLie
The `PauLie` package is an open-source library for studying various algebraic properties of quantum circuits.  This
first release focuses on the classification of the circuit's dynamical Lie algebra given the generators as Paulistrings.

<p align="center">
  <a href="https://qpaulie.github.io/PauLie/">
  <img width=30% src="https://img.shields.io/badge/documentation-blue?style=for-the-badge&logo=read%20the%20docs" alt="Documentation" />
  </a>
</p>


[![Unitary Fund](https://img.shields.io/badge/Supported%20By-UNITARY%20FUND-brightgreen.svg?style=for-the-badge)](https://unitary.fund)

Make sure to have Python >= 3.12 installed.

### Installation
Clone the repository. Once you have [`uv` installed](https://docs.astral.sh/uv/getting-started/installation/#pypi), run:

```sh
uv sync --all-extras --dev
```

### Getting started 

The following code gives an example of usage:

```python
    from paulie.common.pauli_string_factory import get_pauli_string as p
    
    n_qubits = 6
    generators = p(["XYZX", "ZZYZY"], n = n_qubits)
    algebra = generators.get_algebra()
    print(f"number of qubits = {n_qubits}, algebra = {algebra}")
```

Feel free to contribute and check out our open issues. We are also happy to chat with you via [Discord](https://discord.gg/unitary-fund-764231928676089909)

### Testing

In order to run the test suite, run:

```sh
uv run python -m pytest
```

### Type Checking

Check for type errors, improve readability, and assist IDEs, run:

```sh
uv run mypy src/paulie/common src/paulie/classifier/classification.py
```

### Migration from Poetry

This project previously used [Poetry](https://python-poetry.org/) for dependency management.  
As of this release, we have migrated to [`uv`](https://docs.astral.sh/uv/getting-started/installation/#pypi).  

**Old Poetry commands** (for reference during migration):
