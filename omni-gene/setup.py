# -------------------------------------------------------------
# File: setup.py
# Purpose: This script is the build script for setuptools. It tells
#          setuptools about your package (such as the name and version)
#          as well as which code files to include.
# -------------------------------------------------------------
from setuptools import setup, find_packages

# --- Read the version number from a dedicated file ---
# This is a best practice to keep the version number in one place.
version = {}
with open("omni_gene/_version.py") as fp:
    exec(fp.read(), version)


# --- Define the setup configuration ---
setup(
    # The name of your package, this is what users will `pip install`.
    name="omni-gene",
    
    # The version of your package, read from _version.py
    version=version["__version__"],
    
    # Your name and email.
    author="Saket Kumar Jha",
    author_email="saket.jha417@gmail.com",
    
    # A short description of your package.
    description="AI-Powered Multi-Omics Integration and Novel Discovery Platform",
    
    # A long description, usually read from the README.md file.
    long_description=open('README.md').read(),
    long_description_content_type="text/markdown",
    
    # The URL for your project's homepage (e.g., a GitHub repository).
    url="<your_github_repo_url>", 
    
    # Automatically find all packages (directories with an __init__.py file).
    packages=find_packages(),
    
    # Classifiers help users find your project by categorizing it.
    classifiers=[
        "Programming Language :: Python :: 3.9",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    
    # The minimum version of Python required to run your package.
    python_requires='>=3.9',
    
    # A list of other packages that your package depends on.
    # These will be automatically installed when a user installs your package.
    install_requires=[
        "pandas",
        "numpy",
        "scikit-learn",
        "matplotlib",
        "seaborn",
        "scanpy",
        "pysam",
        "PyYAML",
        "rich", # For beautiful terminal output
    ],
    
    # This section creates the command-line tool.
    # It maps the command 'omnigene' to the 'main' function in 'omni_gene/cli.py'.
    entry_points={
        'console_scripts': [
            'omnigene=omni_gene.cli:main',
        ],
    },
)
