import os
import yaml
from datetime import datetime
from ..utils.logger import setup_logger
from .._version import __version__

def run(config):
    """
    Initializes the Omni-Gene pipeline run.
    - Creates a unique, timestamped output directory to prevent overwriting results.
    - Sets up all necessary subdirectories (data, results, images, reports).
    - Configures a logger for this specific run.
    - Returns a comprehensive parameters dictionary for use by other modules.
    """
    params = config.copy()
    
    # Create a unique output directory using the analysis name and current timestamp
    run_timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    analysis_name = params.get('analysis_name', 'omnigene_run')
    run_dir_name = f"{analysis_name}_{run_timestamp}"
    run_dir_path = os.path.join(params['output_dir'], run_dir_name)
    params['run_dir'] = run_dir_path
    
    # Create all necessary subdirectories as specified in the project prompt
    for subdir in ['data', 'results', 'images', 'reports']:
        os.makedirs(os.path.join(run_dir_path, subdir), exist_ok=True)
        
    logger = setup_logger(run_dir_path, analysis_name)
    params['logger'] = logger
    params['version'] = __version__

    logger.info(f"Run directory created at: {run_dir_path}")
    logger.info("Initialization complete.")
    return params
