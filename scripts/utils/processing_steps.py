import os
import scanpy as sc
import pandas as pd
import numpy as np
import scipy
from scipy.sparse import csr_matrix
from itertools import product
from sklearn.neighbors import NearestNeighbors
import time

def write_output(adata, output_adata_file):
    # had to do this because can't write to directory with two different files?
    while True:
        try:
            print(f"Trying to write {output_adata_file}")
            adata.write_h5ad(output_adata_file)
            break
        except Exception as e:
            print(f"An error occurred: {e}")
            time.sleep(5)
        


def hvg(input_adata_file, output_adata_file, hvg_method, num_hvg):
    print(f"input_adata_file: {input_adata_file}, output_adata_file: {output_adata_file}, hvg_method: {hvg_method}, num_hvg: {num_hvg}")
    
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
        print("Finished hvg'")
    write_output(adata, output_adata_file)


def norm(input_adata_file, output_adata_file, norm_method):
    print(f"input_adata_file: {input_adata_file}")
    print(f"output_adata_file: {output_adata_file}")
    print(f"norm_method: {norm_method}")

    # Rest of your function code goes here

    # Check if the output file already exists
    if os.path.exists(output_adata_file):
        print(f"File {output_adata_file} already exists. Skipping normalization.")
        return
        
    adata = sc.read_h5ad(input_adata_file)
    # Apply normalization
    if norm_method == 'norm_log_zscore':
        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)
        sc.pp.scale(adata, max_value=10)
    elif norm_method == 'pearson_residuals':
        sc.experimental.pp.normalize_pearson_residuals(adata)
    elif norm_method == 'pearson_residuals_1':
        sc.experimental.pp.normalize_pearson_residuals(adata, theta=1)
    else:
        raise ValueError("Method not found")
    
    write_output(adata, output_adata_file)
    print("Finished norm")



def hvg_norm(input_adata_file, output_adata_file, hvg_norm_combo, num_hvg):
    print(f"input_adata_file: {input_adata_file}, output_adata_file: {output_adata_file}, combo: {hvg_norm_combo}, num_hvg: {num_hvg}")
    
    # Check if the output file already exists
    if os.path.exists(output_adata_file):
        print(f"File {output_adata_file} already exists. Skipping HVG and Norm.")
        return

    # Load the data
    adata = sc.read_h5ad(input_adata_file)
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=adata.shape[0] * 0.0005) # adata.shape[0] gives the number of cells
    adata_original = adata.copy()

    # Apply HVG selection & normalization method
    if hvg_norm_combo == 'Pearson Residual + Pearson Residual':
        sc.experimental.pp.highly_variable_genes(adata, n_top_genes=num_hvg, flavor='pearson_residuals')
        genes_to_keep = adata.var.highly_variable
        adata = adata_original
        adata = adata[:, genes_to_keep]
        print("Finished HVG with Pearson Residuals")
        sc.experimental.pp.normalize_pearson_residuals(adata)
        print("Finished normalization with Pearson Residuals")
    elif hvg_norm_combo == 'Pearson Residual + norm_log_zscore':
        sc.experimental.pp.highly_variable_genes(adata, n_top_genes=num_hvg, flavor='pearson_residuals')
        genes_to_keep = adata.var.highly_variable
        adata = adata_original
        adata = adata[:, genes_to_keep]
        print("Finished HVG with Pearson Residuals")
        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)
        sc.pp.scale(adata, max_value=10)
        print("Finished normalization with norm_log_zscore")
    elif hvg_norm_combo == 'seurat + norm_log_zscore':
        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, n_top_genes=num_hvg, flavor='seurat')
        genes_to_keep = adata.var.highly_variable
        adata = adata_original
        adata = adata[:, genes_to_keep]
        print("Finished HVG with seurat")
        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)
        sc.pp.scale(adata, max_value=10)
        print("Finished normalization with norm_log_zscore")
    elif hvg_norm_combo == 'seurat + log':
        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, n_top_genes=num_hvg, flavor='seurat')
        genes_to_keep = adata.var.highly_variable
        adata = adata_original
        adata = adata[:, genes_to_keep]
        print("Finished HVG with seurat")
        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)
        print("Finished normalization with log")
    else:
        raise ValueError("Unsupported combo method")

    # Save the processed data
    write_output(adata, output_adata_file)
    print("Output file saved:", output_adata_file)



def pca(input_adata_file, output_adata_file, max_pcs):
    print(f"input_adata_file: {input_adata_file}")
    print(f"output_adata_file: {output_adata_file}")
    print(f"max_pcs: {max_pcs}")
    # Check if the output file already exists
    if os.path.exists(output_adata_file):
        print(f"File {output_adata_file} already exists. Skipping PCA.")
        return
    adata = sc.read_h5ad(input_adata_file)
    sc.tl.pca(adata, n_comps=int(max_pcs))
    adata.write_h5ad(output_adata_file)
    print("Finished pca")

def evaluate(input_adata_file, output_file, label_col, num_nn, num_pcs):
    print(f"input_adata_file: {input_adata_file}")
    print(f"output_file: {output_file}")
    print(f"label_col: {label_col}")
    print(f"num_nn: {num_nn}")
    print(f"num_pcs: {num_pcs}")
    # Rest of your function code goes here

    adata = sc.read_h5ad(input_adata_file)


    results_dict_list = []
    max_num_pcs = adata.obsm['X_pca'].shape[1]
    if num_pcs > max_num_pcs:
        raise ValueError("Not enough PCs to subset")
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
    results_df = results_df.groupby(['label']).mean().reset_index().drop(columns=['cell_index'])
    results_df.to_csv(output_file, sep="\t", index=False)
    print(results_df)







