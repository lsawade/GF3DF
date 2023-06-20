# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information
import os

project = 'GF3DF'
copyright = '2023, Lucas Sawade'
author = 'Lucas Sawade'
release = '0.0.1a'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    # [...]
    'sphinx.ext.intersphinx',
    'sphinx.ext.autodoc',
    'sphinxfortran.fortran_domain',
    'sphinxfortran.fortran_autodoc',
    'sphinx.ext.githubpages',
    'sphinx_design',
    'sphinx_togglebutton',
]
templates_path = ['_templates']
exclude_patterns = ["_build", "Thumbs.db",
                    ".DS_Store", "**.ipynb_checkpoints", "build"]



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'furo'
html_static_path = ['_static']


## -- Options for Sphinx-Fortran ---------------------------------------------
# List of possible extensions in the case of a directory listing
fortran_ext = ['f90', 'F90', 'f95', 'F95']

# This variable must be set with file pattern, like "*.f90", or a list of them.
# It is also possible to specify a directory name; in this case, all files than
# have an extension matching those define by the config variable `fortran_ext`
# are used.
fortran_src = [os.path.abspath('../../src/*.f90'),os.path.abspath('../../src/*.F90')]

# Indentation string or length (default 4). If it is an integer,
# indicates the number of spaces.
fortran_indent = 2