import scanpy as sc
import pandas as pd

def run(adata, params):
    """
    Performs discovery analyses on the integrated multi-omics data.

    This module is the core of the AI-powered discovery phase. It takes the
    unified data representation from the integration step and applies
    statistical and machine learning methods to uncover biological insights.

    This includes:
    - Clustering cells to find novel subtypes based on the combined omics profiles.
    - Running differential expression analysis to find genes that are significantly
      dysregulated between experimental conditions.
    - Identifying the top marker genes that uniquely define each discovered cluster.
    """
    logger = params['logger']
    
    # --- Clustering to find novel cell subtypes ---
    logger.info("Running Leiden clustering on integrated data to find cell subtypes...")
    # The resolution parameter can be tuned in the config file to get more or fewer clusters.
    # This is a key parameter for exploratory analysis.
    sc.tl.leiden(adata, key_added='clusters', resolution=params.get('leiden_resolution', 0.8))
    
    # --- Differential Expression/Abundance Analysis ---
    if params.get('metadata_file'):
        logger.info("Loading metadata for differential analysis...")
        metadata = pd.read_csv(params['metadata_file'], index_col=0)
        # Align metadata with the AnnData object's samples to add condition labels
        adata.obs = adata.obs.join(metadata, how='left')
        
        de_attribute = params.get('de_comparison_attribute')
        # Check if the specified comparison column exists and has more than one group
        if de_attribute and de_attribute in adata.obs.columns and adata.obs[de_attribute].nunique() > 1:
            logger.info(f"Running differential expression analysis on '{de_attribute}'...")
            # We use `use_raw=True` to test on the full, non-subsetted gene set for comprehensive results.
            sc.tl.rank_genes_groups(adata, de_attribute, method='t-test', key_added='de_results', use_raw=True)
            logger.info("Differential expression analysis complete.")
        else:
            logger.warning(f"DE comparison attribute '{de_attribute}' not found in metadata or has fewer than 2 groups. Skipping DE analysis.")
            
    # --- AI Novelty Detection: Find top marker genes for each cluster ---
    # This analysis finds the genes that are most uniquely expressed in each
    # cluster, helping to define the cluster's biological identity.
    logger.info("AI Discovery: Identifying top marker genes that define each cluster...")
    sc.tl.rank_genes_groups(adata, 'clusters', method='t-test', key_added='marker_genes', use_raw=True)
    
    logger.info("Discovery analysis complete.")
    return adata
