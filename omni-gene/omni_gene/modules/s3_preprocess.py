import scanpy as sc
import numpy as np

def run(adata, params):
    """
    Applies tailored pre-processing pipelines to each loaded omics layer.

    This module is responsible for:
    - Applying best-practice quality control (QC) to the RNA data to filter
      out low-quality cells and uninformative genes.
    - Normalizing the RNA data to make counts comparable across cells.
    - Identifying highly variable genes to focus on the most biologically
      meaningful features.
    - Performing basic pre-processing on other loaded modalities, such as
      imputing missing values in proteomics data.
    """
    logger = params['logger']

    # --- Pre-process RNA Layer ---
    logger.info("Applying best-practice preprocessing to RNA layer...")
    # 1. Filter cells with very few genes and genes present in too few cells.
    sc.pp.filter_cells(adata, min_genes=params.get('qc_min_genes', 200))
    sc.pp.filter_genes(adata, min_cells=3)
    
    # 2. Calculate mitochondrial gene percentage for QC.
    #    The gene names for mitochondria can be 'MT-' (human) or 'mt-' (mouse).
    adata.var['mt'] = adata.var_names.str.startswith(('MT-', 'mt-'))
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    
    # 3. Filter out cells with high mitochondrial content (a sign of cell stress/death).
    adata = adata[adata.obs.pct_counts_mt < params.get('qc_max_pct_mt', 20), :].copy()
    
    # 4. Normalize counts to a target sum (e.g., 10,000 reads per cell) and log-transform.
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    
    # 5. Identify highly variable genes to use for downstream analysis.
    sc.pp.highly_variable_genes(adata, n_top_genes=3000)
    
    # Store the full raw data, but subset the main analysis matrix to variable genes.
    adata.raw = adata
    adata = adata[:, adata.var.highly_variable].copy()
    logger.info("RNA preprocessing complete.")

    # --- Pre-process Other Modalities (if they exist) ---
    if 'X_prot' in adata.obsm:
        logger.info("Applying preprocessing to Proteomics layer...")
        # A common step for proteomics is to fill missing values (imputation).
        # This uses a simple mean imputation for demonstration purposes.
        prot_data = adata.obsm['X_prot'].copy()
        col_mean = np.nanmean(prot_data, axis=0)
        # Find the indices of NaN (missing) values
        inds = np.where(np.isnan(prot_data))
        # Replace NaNs with the mean of their respective column (protein)
        prot_data[inds] = np.take(col_mean, inds[1])
        adata.obsm['X_prot'] = prot_data
        logger.info("Proteomics missing value imputation complete.")
        
    return adata
