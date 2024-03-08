import os
import scanpy as sc
import pandas as pd
import numpy as np
import scipy
from scipy.sparse import csr_matrix
from itertools import product
from sklearn.neighbors import NearestNeighbors


def hvg(input_adata_file, output_adata_file, hvg_method, num_hvg):
    # Check if the output file already exists
    if os.path.exists(output_adata_file):
        print(f"File {output_adata_file} already exists. Skipping HVG.")
        return
        
    adata = sc.read_h5ad(input_adata_file)
    sc.pp.filter_genes(adata, min_cells=3)

    adata_original = adata.copy()

    if hvg_method in ['seurat', 'cell_ranger', 'seurat_v3', 'pearson_residuals', 'pearson_residuals_1', 'min_cells', 'total_counts']:
        if hvg_method == 'seurat':
            sc.pp.normalize_total(adata)
            sc.pp.log1p(adata)
            sc.pp.highly_variable_genes(adata, n_top_genes=num_hvg, flavor='seurat')
            genes_to_keep = adata.var.highly_variable
        elif hvg_method == 'cell_ranger':
            sc.pp.normalize_total(adata)
            sc.pp.log1p(adata)
            sc.pp.highly_variable_genes(adata, n_top_genes=num_hvg, flavor='cell_ranger')
            genes_to_keep = adata.var.highly_variable
        elif hvg_method == 'seurat_v3':
            sc.pp.highly_variable_genes(adata, n_top_genes=num_hvg, flavor='seurat_v3')
            genes_to_keep = adata.var.highly_variable
        elif hvg_method == 'pearson_residuals':
            sc.experimental.pp.highly_variable_genes(adata, n_top_genes=num_hvg, flavor='pearson_residuals')
            genes_to_keep = adata.var.highly_variable
        elif hvg_method == 'pearson_residuals_1':
            sc.experimental.pp.highly_variable_genes(adata, n_top_genes=num_hvg, flavor='pearson_residuals', theta=1)
            genes_to_keep = adata.var.highly_variable
            sc.experimental.pp.highly_variable_genes(adata, n_top_genes=num_hvg, flavor='pearson_residuals', theta=1 if hvg_method == 'pearson_residuals_1' else None)
            genes_to_keep = adata.var.highly_variable
        elif hvg_method == 'min_cells':
            sorted_genes = np.argsort(np.sum(adata.X > 0, axis=0))
            genes_to_keep = sorted_genes[-num_hvg:]
        elif hvg_method == 'total_counts':
            sc.pp.normalize_total(adata)
            sorted_genes = np.argsort(np.sum(adata.X, axis=0))
            genes_to_keep = sorted_genes[-num_hvg:]
        else:
            raise ValueError("Method not found")

        adata = adata_original
        
        adata = adata[:, genes_to_keep]
        adata.write_h5ad(output_adata_file)
        print("Finished hvg'")


def norm(input_adata_file, output_adata_file, norm_method):
    # Check if the output file already exists
    if os.path.exists(output_adata_file):
        print(f"File {output_adata_file} already exists. Skipping normalization.")
        return
        
    adata = sc.read_h5ad(input_adata_file)
    # Apply normalization
    if norm_method == 'log_zscore':
        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)
        sc.pp.scale(adata, max_value=10)
    elif norm_method == 'pearson_residuals':
        sc.experimental.pp.normalize_pearson_residuals(adata)
    elif norm_method == 'pearson_residuals_1':
        sc.experimental.pp.normalize_pearson_residuals(adata, theta=1)
    else:
        raise ValueError("Method not found")
    
    adata.write_h5ad(output_adata_file)
    print("Finished norm")



def pca(input_adata_file, output_adata_file, max_pcs):
    # Check if the output file already exists
    if os.path.exists(output_adata_file):
        print(f"File {output_adata_file} already exists. Skipping PCA.")
        return
    adata = sc.read_h5ad(input_adata_file)
    sc.tl.pca(adata, n_comps=int(max_pcs))
    adata.write_h5ad(output_adata_file)
    print("Finished pca")

def evaluate(input_adata_file, output_file, label_col, num_nn, num_pcs_list):
    adata = sc.read_h5ad(input_adata_file)

    if not isinstance(num_pcs_list, list):
        num_pcs_list = [num_pcs_list]  # Ensure it's a list

    results_dict_list = []
    for num_pcs in num_pcs_list:
        X_pca = adata.obsm['X_pca'][:, :num_pcs]
        nbrs = NearestNeighbors(n_neighbors=num_nn, algorithm='brute', n_jobs=-1).fit(X_pca)
        _, knn_indices = nbrs.kneighbors(X_pca)


        for label_value in adata.obs[label_col].unique():
            cells_with_label_idx = np.where(adata.obs[label_col] == label_value)[0]

            for i in cells_with_label_idx:
                neighbors = knn_indices[i]
                neighbors_count = np.sum(adata.obs[label_col][neighbors] == label_value)

                result_dict = {
                    'num_PCs': num_pcs,
                    'label': label_value,
                    'cell_index': i,
                    'neighbors_count': neighbors_count,
                }
                results_dict_list.append(result_dict)

    results_df = pd.DataFrame(results_dict_list)
    results_df = results_df.groupby(['num_PCs', 'label']).mean().reset_index().drop(columns=['cell_index'])
    results_df.to_csv(output_file, sep="\t", index=False)







