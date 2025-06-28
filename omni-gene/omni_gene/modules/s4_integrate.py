import scanpy as sc
import numpy as np

def run(adata, params):
    """
    Integrates multi-omics data using the concatenated PCA method.
    
    This approach creates a unified feature space that represents all
    omics layers simultaneously. It works by:
    1. Performing Principal Component Analysis (PCA) on each modality separately
       to reduce its dimensionality while retaining its main sources of variation.
    2. Concatenating (stacking) these individual PCA results into a single,
       wide matrix. This new matrix represents the integrated state of the cell.
    3. Running downstream steps like UMAP and clustering on this integrated representation.
    """
    logger = params['logger']
    pca_results = []
    
    # --- RNA PCA ---
    logger.info("Running PCA on RNA data...")
    # Scale the data (z-score) before PCA. Required for linear dimensionality reduction.
    sc.pp.scale(adata)
    sc.tl.pca(adata, n_comps=30, svd_solver='arpack')
    pca_results.append(adata.obsm['X_pca'])
    
    # --- PCA on Other Modalities ---
    # Loop through other potential data types that might have been loaded.
    for mod_name in ['prot', 'geno']:
        if f'X_{mod_name}' in adata.obsm:
            logger.info(f"Running PCA on {mod_name.capitalize()} data...")
            # Create a temporary AnnData object just for running the PCA algorithm.
            temp_adata = sc.AnnData(adata.obsm[f'X_{mod_name}'])
            sc.pp.scale(temp_adata)
            sc.tl.pca(temp_adata, n_comps=15, svd_solver='arpack')
            pca_results.append(temp_adata.obsm['X_pca'])
            
    # --- Concatenate PCA results ---
    logger.info("Concatenating PCA embeddings into a single multi-modal representation.")
    # np.hstack stacks arrays in sequence horizontally (column-wise).
    adata.obsm['X_integrated'] = np.hstack(pca_results)
    
    # --- Compute UMAP on the integrated embedding ---
    logger.info("Computing neighbors and UMAP on integrated data...")
    # This step is crucial for visualization and clustering on the combined data.
    # We explicitly tell scanpy to use our new 'X_integrated' representation.
    sc.pp.neighbors(adata, use_rep='X_integrated')
    sc.tl.umap(adata)
    
    logger.info("Multi-omics integration complete.")
    return adata
