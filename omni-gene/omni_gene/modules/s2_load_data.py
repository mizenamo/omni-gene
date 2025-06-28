import os
import scanpy as sc
import pandas as pd
import numpy as np
import allel  # The correct library for VCF handling

def run(params):
    """
    Loads, validates, and harmonizes all specified multi-omics data types.

    This rewritten version includes:
    - More robust VCF file parsing using scikit-allel.
    - Improved logic and error reporting for harmonizing samples across all modalities.
    """
    logger = params['logger']
    modalities = {}

    # --- Step 1: Load each modality individually ---

    # Load RNA data
    if 'RNA' in params['data_paths']:
        rna_path = params['data_paths']['RNA']
        logger.info(f"Loading Transcriptome data from: {rna_path}")
        try:
            if os.path.isdir(rna_path):
                modalities['rna'] = sc.read_10x_mtx(rna_path, var_names='gene_symbols', cache=True)
            elif rna_path.endswith('.h5ad'):
                modalities['rna'] = sc.read_h5ad(rna_path)
            elif rna_path.endswith('.h5'):
                modalities['rna'] = sc.read_10x_h5(rna_path)
            else:
                logger.error(f"Unsupported RNA data format or path: {rna_path}")
            if 'rna' in modalities:
                logger.info(f"-> Found RNA data with {modalities['rna'].n_obs} cells/samples. Sample IDs begin with: {modalities['rna'].obs.index[:3].tolist()}...")
        except Exception as e:
            logger.error(f"Failed to load RNA data. Error: {e}")

    # Load Proteomics data
    if 'Proteomics' in params['data_paths']:
        prot_path = params['data_paths']['Proteomics']
        logger.info(f"Loading Proteome data from: {prot_path}")
        try:
            prot_df = pd.read_csv(prot_path, index_col=0)
            adata_prot = sc.AnnData(prot_df.T)
            modalities['prot'] = adata_prot
            logger.info(f"-> Found Proteomics data with {modalities['prot'].n_obs} samples. Sample IDs are: {modalities['prot'].obs.index[:3].tolist()}...")
        except Exception as e:
            logger.error(f"Failed to load Proteomics data. Error: {e}")
            
    # Load Genomics data
    if 'Genomics' in params['data_paths']:
        vcf_path = params['data_paths']['Genomics']
        logger.info(f"Loading Genome data from: {vcf_path}")
        try:
            vcf_data = allel.read_vcf(vcf_path, fields=['samples', 'calldata/GT'])
            if vcf_data and isinstance(vcf_data['calldata/GT'], allel.GenotypeArray):
                genotype_matrix = vcf_data['calldata/GT'].count_alleles().to_n_alt().T
                adata_geno = sc.AnnData(genotype_matrix)
                adata_geno.obs_names = vcf_data['samples']
                modalities['geno'] = adata_geno
                logger.info(f"-> Found Genomics data with {modalities['geno'].n_obs} samples. Sample IDs begin with: {modalities['geno'].obs.index[:3].tolist()}...")
            else:
                logger.warning("VCF file did not contain expected genotype data ('calldata/GT'). Skipping Genomics.")
        except Exception as e:
            logger.error(f"Failed to load or process VCF file. Error: {e}")

    # --- Step 2: Harmonize Samples (The CRITICAL FIX) ---

    if not modalities or 'rna' not in modalities:
        raise ValueError("Pipeline failed: RNA data is required to serve as the base for integration.")
        
    logger.info("Harmonizing samples across all loaded modalities...")
    
    # Start with the sample names from the base RNA modality
    common_samples = modalities['rna'].obs.index
    
    # Iteratively find the intersection of sample names across all other modalities
    for mod_name, mod_data in modalities.items():
        if mod_name != 'rna':
            common_samples = common_samples.intersection(mod_data.obs.index)

    if len(common_samples) == 0:
        # Improved Error Reporting
        rna_samples = modalities.get('rna', sc.AnnData()).obs.index.tolist()
        prot_samples = modalities.get('prot', sc.AnnData()).obs.index.tolist()
        geno_samples = modalities.get('geno', sc.AnnData()).obs.index.tolist()
        error_msg = (
            "Pipeline failed: No common samples found between the loaded omics datasets.\n"
            f"Please check that the sample IDs match between your files.\n"
            f"  - RNA samples found: {rna_samples[:5]}...\n"
            f"  - Proteomics samples found: {prot_samples[:5]}...\n"
            f"  - Genomics samples found: {geno_samples[:5]}..."
        )
        raise ValueError(error_msg)

    logger.info(f"Found {len(common_samples)} common samples across all datasets. Subsetting all data to these samples.")

    # --- Step 3: Create the final, harmonized AnnData object ---
    
    # Filter the base RNA object first
    adata = modalities['rna'][common_samples, :].copy()
    adata.obs_names_make_unique()

    # Filter other modalities and add them to the .obsm slot
    for mod_name, mod_data in modalities.items():
        if mod_name != 'rna':
            aligned_mod_data = mod_data[common_samples, :].copy()
            adata.obsm[f'X_{mod_name}'] = aligned_mod_data.X

    logger.info(f"All data loaded and harmonized for {adata.n_obs} common samples.")
    return adata
