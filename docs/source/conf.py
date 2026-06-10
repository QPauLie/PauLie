"""
Configuration file for the Sphinx documentation builder.
"""
# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import os
import sys
sys.path.insert(0, os.path.abspath('../../'))

import matplotlib  # pylint: disable=wrong-import-position
matplotlib.use("Agg")  # Render animations headless during the documentation build.

from intersphinx_registry import get_intersphinx_mapping  # pylint: disable=wrong-import-position

project = 'paulie'
copyright = '2026, PauLie contributors'  # pylint: disable=redefined-builtin
author = 'PauLie contributors'
release = '0.0.1'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
'sphinx.ext.autodoc',
'sphinx.ext.autosummary',
'sphinx.ext.napoleon',
'sphinx.ext.mathjax',
'sphinx.ext.intersphinx',
'sphinxcontrib.bibtex',
]

bibtex_bibfiles = ["references.bib"]
bibtex_reference_style = "author_year"
autosummary_generate = True
autosummary_generate_overwrite = True
autodoc_typehints = "none"
napoleon_google_docstring = True


templates_path = ['_templates']
exclude_patterns = []

# -----------------------------------------------------------------------------
# Intersphinx configuration
# -----------------------------------------------------------------------------
intersphinx_mapping = get_intersphinx_mapping(
    packages={"python", "numpy"}
)

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "pydata_sphinx_theme"
html_title = "PauLie"

html_static_path = ['_static']
html_css_files = ['custom.css']
mathjax_path = 'https://cdn.jsdelivr.net/npm/mathjax@4/tex-mml-chtml.js'
html_extra_path=['media']


# -- Classification animations -----------------------------------------------
# The interactive animations embedded in user/classification.rst are generated at build time
# into the media directory rather than committed to the repository.

def _generate_animations(_app) -> None:
    """
    Generate the embedded classification animations into the media directory.

    Args:
        _app (sphinx.application.Sphinx): Sphinx application (unused).
    Returns:
        None
    """
    from paulie import (  # pylint: disable=import-outside-toplevel
        animation_anti_commutation_graph,
        get_pauli_string as p,
    )
    media_dir = os.path.join(os.path.dirname(__file__), "media")
    os.makedirs(media_dir, exist_ok=True)
    examples = {
        "example_c": p(["IYZI", "IIXX", "IIYZ", "IXXI", "XXII", "YZII"]),
        "example_d": p(["XY", "XZ"], n=4),
    }
    for name, generators in examples.items():
        anim = animation_anti_commutation_graph(generators, interval=1200)
        html = anim.to_jshtml(default_mode="loop")
        with open(os.path.join(media_dir, f"{name}.html"), "w", encoding="utf-8") as file:
            file.write(html)


def setup(app) -> None:
    """
    Connect the animation generation to the documentation build.

    Args:
        app (sphinx.application.Sphinx): Sphinx application.
    Returns:
        None
    """
    app.connect("builder-inited", _generate_animations)
