# PauLie
The `PauLie` package is an open-source library for studying various algebraic properties of quantum circuits.  This
first release focuses on the classification of the circuit's dynamical Lie algebra given the generators as Paulistrings.

[![Unitary Fund](https://img.shields.io/badge/Supported%20By-UNITARY%20FUND-brightgreen.svg?style=for-the-badge)](https://unitary.fund)

Make sure to have Python >= 3.12 installed.

### Installation
Once you have `poetry` installed, run:

```sh
poetry install
poetry shell
```
Notes for Windows.
If there are errors installing dependent packages, then you need to configure the poetry installation environment in the folder.
```sh
poetry config virtualenvs.in-project true
```
There is no need to launch poetry shell.
To activate use
```sh
poetry env activate
```

### Description

#### 1. src
Package source code

#### 2. docs
sphinx documentation
https://qpaulie.github.io/PauLie/

#### 3. examples
Examples of use

#### 4. tests
Tests for pytest

Feel free to contribute and check out our open issues. 





