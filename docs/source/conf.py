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

project = 'paulie'
copyright = '2024, PauLie contributors'
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
]

autosummary_generate = True
autosummary_generate_overwrite = True
autodoc_typehints = "none"
napoleon_google_docstring = True
napoleon_numpy_docstring = True

templates_path = ['_templates']
exclude_patterns = []



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "pydata_sphinx_theme"
html_title = "PauLie"

html_static_path = ['_static']
html_css_files = ['custom.css']
mathjax_path="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML"
html_extra_path=['media']
