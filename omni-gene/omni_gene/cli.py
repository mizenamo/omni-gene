import argparse
from .main_pipeline import run_omni_gene_pipeline
from ._version import __version__

def main():
    """
    Defines and parses the command-line interface for the Omni-Gene tool.

    This function sets up the user-facing command, 'omnigene', and defines
    the arguments it accepts, such as the required configuration file. It then
    calls the main pipeline orchestrator to start the analysis.
    """
    parser = argparse.ArgumentParser(
        description=f"Omni-Gene v{__version__}: AI-Powered Multi-Omics Integration and Discovery Platform.",
        # RawTextHelpFormatter allows for better formatting of the help message.
        formatter_class=argparse.RawTextHelpFormatter
    )
    
    # Define the required '--config' argument
    parser.add_argument(
        '-c', '--config',
        type=str,
        required=True,
        help='Path to the main configuration YAML file (e.g., config.yaml).'
    )
    
    # Define an optional '--version' argument to display the tool's version
    parser.add_argument(
        '--version',
        action='version',
        version=f'%(prog)s {__version__}',
        help="Show the program's version number and exit."
    )
    
    # Parse the arguments provided by the user in the terminal
    args = parser.parse_args()
    
    # Launch the main analysis pipeline, passing the path to the config file
    run_omni_gene_pipeline(config_path=args.config)

