# -------------------------------------------------------------
# File: omni_gene/utils/__init__.py
# Purpose: This file makes the 'utils' directory a Python package.
#          By importing key utility functions here, they become
#          easier to access from other parts of the application.
# -------------------------------------------------------------

# Import the setup_logger function from the logger module.
# This allows other scripts to import it directly from the utils package,
# for example: from omni_gene.utils import setup_logger
from .logger import setup_logger