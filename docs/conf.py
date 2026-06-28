# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# http://www.sphinx-doc.org/en/master/config

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
import datetime
from importlib.metadata import version as get_version, PackageNotFoundError

sys.path.insert(0, os.path.abspath("../src"))


# The master toctree document.
master_doc = "index"

# -- Project information -----------------------------------------------------

project = "pyMARS"
author = "Phillip Mestas, Parker Clayton, and Kyle Niemeyer"
this_year = datetime.date.today().year
copyright = "{}, {}".format(this_year, author)

# The version info for the project you're documenting, acts as replacement for
# |version| and |release|, also used in various other places throughout the
# built documents.
#
# The full version, including alpha/beta/rc tags. Looked up by the PyPI
# distribution name ("nrg-pymars"); the import package is still "pymars".
try:
    release = get_version("nrg-pymars")
except PackageNotFoundError:
    release = "unknown"
# The short X.Y version.
version = ".".join(release.split(".")[:2])


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.doctest",
    "sphinx.ext.napoleon",
    "sphinx.ext.mathjax",
    "sphinx.ext.intersphinx",
    "sphinxcontrib.bibtex",
]

autodoc_default_options = {"members": True}
autoclass_content = "class"
napoleon_numpy_docstring = True
napoleon_google_docstring = False

# sphinxcontrib-bibtex (v2+) requires bib files to be declared here, and the
# citation style is set globally rather than per ``bibliography`` directive.
bibtex_bibfiles = ["refs.bib"]
bibtex_default_style = "unsrt"

intersphinx_mapping = {
    "python": ("https://docs.python.org/3.14", None),
    "numpy": ("https://docs.scipy.org/doc/numpy/", None),
    "scipy": ("https://docs.scipy.org/doc/scipy/", None),
    "networkx": ("https://networkx.github.io/documentation/stable/", None),
    "cantera": ("https://cantera.org/stable/", None),
    #'pytables': ('http://www.pytables.org/usersguide/libref/', None)
}

# Make the module index more useful by sorting on the module name
# instead of the package name
modindex_common_prefix = ["pymars."]

# Suppress duplicate citation warnings. "ref.citation" covers the docutils-level
# warnings (and, as a side effect, undefined-citation and unused-reference
# warnings). "bibtex.duplicate_citation" covers a reference cited from more than
# one page-local bibliography (e.g., niemeyer:2015 appears in both theory and
# psr_theory); label collisions are already avoided via the bibliography
# directive's ``:labelprefix:``.
suppress_warnings = ["ref.citation", "bibtex.duplicate_citation"]

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store", "_graph"]

source_suffix = {
    ".rst": "restructuredtext",
    #  '.txt': 'markdown',
    #  '.md': 'markdown',
}

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = "sphinx"

# If true, `todo` and `todoList` produce output, else they produce nothing.
todo_include_todos = True


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "alabaster"

# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
html_theme_options = {
    "github_user": "Niemeyer-Research-Group",
    "github_repo": "pyMARS",
    "github_banner": True,
    "github_button": True,
    "show_powered_by": True,
    "logo": "pymars-logo.png",
}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]

# -- Options for LaTeX output ---------------------------------------------

latex_elements = {
    # The paper size ('letterpaper' or 'a4paper').
    #
    # 'papersize': 'letterpaper',
    # The font size ('10pt', '11pt' or '12pt').
    #
    # 'pointsize': '10pt',
    # Additional stuff for the LaTeX preamble.
    #
    # 'preamble': '',
    # Latex figure (float) alignment
    #
    # 'figure_align': 'htbp',
}

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title,
#  author, documentclass [howto, manual, or own class]).
latex_documents = [
    (master_doc, "pyMARS.tex", "pyMARS Documentation", author, "manual"),
]


# -- Options for manual page output ---------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [(master_doc, "pymars", "pyMARS Documentation", [author], 1)]


# -- Options for Texinfo output -------------------------------------------

# Grouping the document tree into Texinfo files. List of tuples
# (source start file, target name, title, author,
#  dir menu entry, description, category)
texinfo_documents = [
    (
        master_doc,
        "pyMARS",
        "pyMARS Documentation",
        author,
        "pyMARS",
        "One line description of project.",
        "Miscellaneous",
    ),
]
