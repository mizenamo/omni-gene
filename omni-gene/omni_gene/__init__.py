# -------------------------------------------------------------
# File: omni_gene/__init__.py
# Purpose: This file makes the 'omni_gene' directory a Python package.
#          By importing the main function and version here, they can be
#          easily accessed by other developers using your package.
# -------------------------------------------------------------

# Expose the main pipeline function at the package level.
# This allows programmatic use, e.g., `from omni_gene import run_omni_gene_pipeline`.
from .main_pipeline import run_omni_gene_pipeline

# Expose the version number at the package level.
# This allows users to check the installed version, e.g., `import omni_gene; print(omni_gene.__version__)`.
from ._version import __version__
