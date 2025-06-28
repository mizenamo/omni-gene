import scanpy as sc
import pandas as pd
import os
import zipfile
from datetime import datetime

def run(adata, params):
    """
    Generates all final output files as specified in the project prompt
    and creates a single, downloadable ZIP archive of the entire run.

    This module is responsible for:
    - Saving high-resolution UMAP plots for clusters and conditions.
    - Saving key differential expression plots (e.g., volcano, dot plots).
    - Writing the final processed AnnData object to a file for later use.
    - Saving statistical results like DE tables and marker gene lists to CSV files.
    - Creating a summary README file for the analysis run.
    - Compiling all outputs into a single, compressed ZIP file.
    """
    logger = params['logger']
    run_dir = params['run_dir']
    img_dir = os.path.join(run_dir, 'images')
    res_dir = os.path.join(run_dir, 'results')
    data_dir = os.path.join(run_dir, 'data')
    
    # --- Generate and Save Plots ---
    logger.info("Generating final plots for the report...")
    
    # Save UMAP of integrated clusters
    sc.pl.umap(adata, color='clusters', save='_clusters.png', show=False, legend_loc='on data', title='Integrated Cell Subtypes')
    
    # Save UMAP colored by experimental condition, if it exists
    de_attribute = params.get('de_comparison_attribute')
    if de_attribute and de_attribute in adata.obs.columns:
        sc.pl.umap(adata, color=de_attribute, save=f'_{de_attribute}.png', show=False, title=f'Samples by {de_attribute}')
    
    # Save differential expression results plots
    if 'de_results' in adata.uns:
        sc.pl.rank_genes_groups(adata, n_genes=20, sharey=False, key='de_results', save='_de_summary.png', show=False)
        sc.pl.rank_genes_groups_volcano(adata, n_genes=10, key='de_results', save='_volcano.png', show=False)

    # Save marker gene plots
    if 'marker_genes' in adata.uns:
        sc.pl.rank_genes_groups_dotplot(adata, n_genes=5, key='marker_genes', save='_marker_genes.png', show=False)

    # Move all generated plots from the default location to the 'images' directory
    for f in os.listdir('.'):
        if f.startswith(('umap_', 'rank_genes_groups_')) and f.endswith('.png'):
            try:
                os.rename(f, os.path.join(img_dir, f.lstrip('figures/')))
            except FileNotFoundError:
                logger.warning(f"Could not find plot file {f} to move, it may not have been generated.")

    # --- Generate Data and Results Tables ---
    logger.info("Generating results tables...")
    # 1. Save the final, processed AnnData object. This is crucial for reproducibility.
    adata.write(os.path.join(data_dir, 'final_processed_object.h5ad'))
    
    # 2. Save differential expression results table
    if 'de_results' in adata.uns:
        de_df = sc.get.rank_genes_groups_df(adata, group=None, key='de_results')
        de_df.to_csv(os.path.join(res_dir, 'differential_expression_results.csv'), index=False)
        
    # 3. Save marker gene results table
    if 'marker_genes' in adata.uns:
        marker_df = sc.get.rank_genes_groups_df(adata, group=None, key='marker_genes')
        marker_df.to_csv(os.path.join(res_dir, 'cluster_marker_genes.csv'), index=False)
        
    # --- Create Final Report and ZIP Archive ---
    readme_path = os.path.join(run_dir, 'README.txt')
    logger.info("Creating summary README file...")
    with open(readme_path, 'w') as f:
        f.write(f"Omni-Gene Analysis Report\n")
        f.write("============================\n\n")
        f.write(f"Analysis Name: {params['analysis_name']}\n")
        f.write(f"Run Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Tool Version: {params['version']}\n\n")
        f.write("This archive contains all data, results, and images from the pipeline run.\n\n")
        f.write("Directory Structure:\n")
        f.write("- /data: Contains the final processed AnnData object (.h5ad) for re-analysis.\n")
        f.write("- /results: Contains all statistical output tables (e.g., DE genes, marker genes) in CSV format.\n")
        f.write("- /images: Contains all high-resolution plots and graphs in PNG format.\n")

    zip_path = f"{run_dir}.zip"
    logger.info(f"Creating final ZIP archive at: {zip_path}")
    with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
        for root, _, files in os.walk(run_dir):
            for file in files:
                file_path = os.path.join(root, file)
                # Create a relative path for the files inside the zip archive
                archive_name = os.path.relpath(file_path, params['output_dir'])
                zipf.write(file_path, arcname=archive_name)

    logger.info("Report generation complete. Your results are available in the ZIP archive.")
