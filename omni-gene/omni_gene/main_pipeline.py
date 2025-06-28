import yaml
from .modules import (
    s1_initialize,
    s2_load_data,
    s3_preprocess,
    s4_integrate,
    s5_discover,
    s6_report
)

def run_omni_gene_pipeline(config_path):
    """
    Orchestrates the entire Omni-Gene end-to-end workflow, from setup
    to final report generation.
    """
    # Load the user-defined parameters from the YAML config file.
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)

    # --- Step 1: Setup and Initialization ---
    params = s1_initialize.run(config)
    logger = params['logger']
    logger.info(f"Starting Omni-Gene v{params['version']} pipeline.")

    # --- Step 2: Load, Validate, and Harmonize Data ---
    adata = s2_load_data.run(params)

    # --- Step 3: Pre-process Each Omics Layer ---
    adata = s3_preprocess.run(adata, params)
    
    # --- Step 4: Integrate Multi-Omics Data ---
    adata = s4_integrate.run(adata, params)
    
    # --- Step 5: AI-Powered Discovery and Functional Analysis ---
    adata = s5_discover.run(adata, params)

    # --- Step 6: Generate Final Report and Archive ---
    s6_report.run(adata, params)
    
    logger.info("Omni-Gene pipeline finished successfully!")
