import os
import sys
sys.path.insert(0, os.path.abspath('../'))

# -- Project information -----------------------------------------------------

project = 'ChromProcess'
copyright = '2022, W. E. Robinson'
author = 'W. E. Robinson'

# The full version, including alpha/beta/rc tags
release = '0.0.0'

# -- General configuration ---------------------------------------------------

extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'numpydoc']

templates_path = ['_templates']

exclude_patterns = []

autoapi_dirs = ['../ChromProcess']

# -- Options for HTML output -------------------------------------------------

html_theme = 'pyramid'

html_static_path = ['_static']

# -- numpydoc config -------------------------------------------------

numpydoc_use_plots = False

numpydoc_show_class_members = True

numpydoc_show_inherited_class_members = False

numpydoc_class_members_toctree = False

numpydoc_citation_re = '[\w-]+'

numpydoc_use_blockquotes = True

numpydoc_attributes_as_param_list = True

numpydoc_xref_param_type = True

#numpydoc_xref_aliases
#numpydoc_xref_ignore
#numpydoc_validation_checks
#numpydoc_validation_excluden
