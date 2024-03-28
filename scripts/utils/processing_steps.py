import os
import scanpy as sc
import pandas as pd
import numpy as np
import scipy
from scipy.sparse import csr_matrix
from itertools import product
from sklearn.neighbors import NearestNeighbors
import time
import gzip
import shutil

def write_output(adata, output_file):
    # had to do this because can't write to directory with two different files?
    while True:
        try:
            print(f"Trying to write {output_file}")
            adata.write_h5ad(output_file)
            break
        except Exception as e:
            print(f"An error occurred: {e}")
            time.sleep(5)
        

def clean(input_adata_file, output_file, max_num_hvg_out_file):
    """Prepare 'cleaned_input.h5ad' from the input AnnData object."""
    adata = sc.read_h5ad(input_adata_file)
    adata.var_names_make_unique()
    if not isinstance(adata.X, scipy.sparse.spmatrix):
        adata.X = csr_matrix(adata.X)
    adata.layers['raw'] = adata.X
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=5)
    adata.write_h5ad(output_file)
    print("Finished preparing 'cleaned_input.h5ad'")

    with open(max_num_hvg_out_file, 'w') as f:
        f.write(str(adata.shape[1]))


def clean_CITE_seq(input_adata_file, output_file, protein_processed_out_file, max_num_hvg_out_file):
    adata = sc.read_h5ad(input_adata_file)
    adata.var_names_make_unique()
    rna_adata = adata[:, adata.var['feature_types'] == 'Gene Expression']
    protein_adata = adata[:, adata.var['feature_types'] == 'Antibody Capture']
    if not isinstance(rna_adata.X, scipy.sparse.spmatrix):
        rna_adata.X = csr_matrix(rna_adata.X)
    rna_adata.layers['raw'] = rna_adata.X
    sc.pp.filter_cells(rna_adata, min_genes=200)
    sc.pp.filter_genes(rna_adata, min_cells=5)
    rna_adata.write_h5ad(output_file)

    print("Finished preparing 'cleaned_input.h5ad'")
    cells_in_rna_adata = rna_adata.obs.index
    protein_adata = protein_adata[cells_in_rna_adata, :]
    clr_normalize(protein_adata)
    sc.pp.scale(protein_adata, max_value=10)
    protein_adata.write_h5ad(protein_processed_out_file)

    with open(max_num_hvg_out_file, 'w') as f:
        f.write(str(rna_adata.shape[1]))


def hvg(input_adata_file, output_file, hvg_method, num_hvg):
    print(f"input_adata_file: {input_adata_file}, output_file: {output_file}, hvg_method: {hvg_method}, num_hvg: {num_hvg}")
    
    # Check if the output file already exists
    if os.path.exists(output_file):
        print(f"File {output_file} already exists. Skipping HVG.")
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
    write_output(adata, output_file)


def norm(input_adata_file, output_file, norm_method):
    print(f"input_adata_file: {input_adata_file}")
    print(f"output_file: {output_file}")
    print(f"norm_method: {norm_method}")

    # Rest of your function code goes here

    # Check if the output file already exists
    if os.path.exists(output_file):
        print(f"File {output_file} already exists. Skipping normalization.")
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
    
    write_output(adata, output_file)
    print("Finished norm")



def hvg_norm(input_adata_file, output_file, hvg_norm_combo, num_hvg):
    print(f"input_adata_file: {input_adata_file}, output_file: {output_file}, combo: {hvg_norm_combo}, num_hvg: {num_hvg}")
    
    # Check if the output file already exists
    if os.path.exists(output_file):
        print(f"File {output_file} already exists. Skipping HVG and Norm.")
        return

    # Load the data
    adata = sc.read_h5ad(input_adata_file)
    adata_original = adata.copy()

    # Convert num_hvg to integer if it is a string
    if isinstance(num_hvg, str):
        if num_hvg.isdigit():
            num_hvg = int(num_hvg)
        elif num_hvg == 'max_num_hvg':
            num_hvg = adata.shape[1]
        else:
            raise ValueError("num_hvg must be an integer or 'max_num_hvg'")

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
    elif hvg_norm_combo == 'seurat + norm_log':
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
    elif hvg_norm_combo == 'sctransform':
        print("Running sctransform")
        print("Number of cells, gene", adata.shape)
        # Example usage
        input_dir = output_file + "_sctransform_input"
        output_dir = output_file + "_sctransform_output"
        write_10x_mtx_gz(input_dir, adata)
        cmd = f"""
        /Genomics/pritykinlab/dillon/software/miniconda/envs/envs/seurat/bin/Rscript utils/sctransform.R {input_dir} {output_dir} {num_hvg}
        """
        if os.path.exists(output_dir):
            # remove the output directory if it already exists
            os.system(f"rm -r {output_dir}")
        exit_status = os.system(cmd)

        if exit_status != 0:
            raise Exception("R script execution was not successful.")

        time.sleep(60) # wait for everything to be written
        adata = read_10x_mtx_gz(output_dir)
        adata.obs = adata_original.obs
    else:
        raise ValueError("Unsupported combo method")

    # Save the processed data
    write_output(adata, output_file)
    print("Output file saved:", output_file)

def write_10x_mtx_gz(output_dir, adata):
    """
    Writes the AnnData object to files compatible with the 10x format and compresses them with gzip.
    Parameters:
    - output_dir: Directory to write the 10x files to.
    - adata: The AnnData object.
    """
    print("Writing 10x File")
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Extract the gene expression matrix (in COO format) and transpose it
    mat = adata.X.T

    # Convert the matrix to integer type if it's not already
    if mat.dtype != 'int':
        mat = mat.astype('int')

    # Check if the matrix is sparse, and convert to COO format if it's not
    if not isinstance(mat, scipy.sparse.coo_matrix):
        mat = scipy.sparse.coo_matrix(mat)

    # Write the matrix to a .mtx file and then compress it
    mtx_path = os.path.join(output_dir, 'matrix.mtx')
    scipy.io.mmwrite(mtx_path, mat)
    with open(mtx_path, 'rb') as f_in:
        with gzip.open(mtx_path + '.gz', 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    os.remove(mtx_path)

    # Function to compress a file
    def compress_file(input_path):
        with open(input_path, 'rb') as f_in:
            with gzip.open(input_path + '.gz', 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        os.remove(input_path)

    # Write genes.tsv and compress it
    genes = adata.var_names
    cleaned_genes = []

    gene_counts = {}
    for gene in adata.var_names:
        # Replace newline with double dashes and underscore with dash
        new_gene = gene.replace('\n', '--').replace('_', '-')

        # Count occurrences and append a number if there's a duplicate
        if new_gene in gene_counts:
            gene_counts[new_gene] += 1
            new_gene += f".{gene_counts[new_gene]}"
        else:
            gene_counts[new_gene] = 0

        cleaned_genes.append(new_gene)

    genes = pd.DataFrame({'gene_ids': cleaned_genes})
    genes_path = os.path.join(output_dir, 'features.tsv')
    genes['genes'] = genes['gene_ids']
    genes.to_csv(genes_path, sep='\t', index=False, header=False)
    compress_file(genes_path)

    # Write barcodes.tsv and compress it
    barcodes = pd.DataFrame({'barcodes': adata.obs_names})
    barcodes_path = os.path.join(output_dir, 'barcodes.tsv')
    barcodes.to_csv(barcodes_path, sep='\t', index=False, header=False)
    compress_file(barcodes_path)
    print("Finisehd Writing 10x files (this is in python)")

from anndata import AnnData

def read_10x_mtx_gz(input_dir):
    # Read the matrix
    import h5py

    with h5py.File(os.path.join(input_dir, 'matrix.h5'), 'r') as file:
        matrix = file['count_matrix'][:]

    genes = pd.read_csv(os.path.join(input_dir, 'features.tsv.gz'), names=['gene_symbol'], header=None, sep='\t',index_col=0)
    barcodes = pd.read_csv(os.path.join(input_dir, 'barcodes.tsv.gz'), names=['cell_barcode'], header=None, index_col=0)

    # Create the AnnData object
    adata = AnnData(X=matrix.T, obs=barcodes, var=genes)
    print(adata)
    return adata



def pca(input_adata_file, output_file, max_pcs):
    print(f"input_adata_file: {input_adata_file}")
    print(f"output_file: {output_file}")
    print(f"max_pcs: {max_pcs}")
    # Check if the output file already exists
    if os.path.exists(output_file):
        print(f"File {output_file} already exists. Skipping PCA.")
        return
    adata = sc.read_h5ad(input_adata_file)
    sc.tl.pca(adata, n_comps=int(max_pcs))
    adata.write_h5ad(output_file)
    print("Finished pca")

def evaluate_old(input_adata_file, output_file, label_col, num_nn, num_pcs):
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
                'label': label_value,
                'cell_index': i,
                'neighbors_count': neighbors_count,
            }
            results_dict_list.append(result_dict)

    results_df = pd.DataFrame(results_dict_list)
    results_df = results_df.groupby(['label']).mean().reset_index().drop(columns=['cell_index'])
    results_df.to_csv(output_file, sep="\t", index=False)

def evaluate(input_adata_file, output_file, label_col, num_nn, num_pcs):
    print(f"input_adata_file: {input_adata_file}")
    print(f"output_file: {output_file}")
    print(f"label_col: {label_col}")
    print(f"num_nn: {num_nn}")
    print(f"num_pcs: {num_pcs}")

    adata = sc.read_h5ad(input_adata_file)

    # Precompute the label counts
    label_counts = adata.obs[label_col].value_counts()
    sufficient_labels = label_counts > num_nn

    # Precompute indices for each label
    label_indices = {label: (adata.obs[label_col] == label) for label in label_counts.index}

    # Check for sufficient PCs
    max_num_pcs = adata.obsm['X_pca'].shape[1]
    if num_pcs > max_num_pcs:
        raise ValueError("Not enough PCs to subset")
    X_pca = adata.obsm['X_pca'][:, :num_pcs]


    results_dict_list = []

    for cell_idx in range(adata.shape[0]):
        label = adata.obs[label_col][cell_idx]
        # if label is nan just skp
        if pd.isna(label):
            continue

        # Skip labels that do not have sufficient instances
        if not sufficient_labels[label]:
            continue

        same_label_indices = label_indices[label]
        diff_label_indices = ~same_label_indices

        distances = np.linalg.norm(X_pca - X_pca[cell_idx], axis=1)
        same_label_dists = distances[same_label_indices]
        diff_label_dists = distances[diff_label_indices]

        nth_neighbor_dist = same_label_dists[num_nn]

        num_closer_diff_label = np.sum(diff_label_dists < nth_neighbor_dist)

        num_cells_of_label = label_counts[label]
        num_cells_of_diff_label = adata.shape[0] - num_cells_of_label
        enrichment = num_closer_diff_label / num_cells_of_diff_label * num_cells_of_label
        enrichment = enrichment / num_nn

        results_dict_list.append({
            'label': label,
            'cell_index': cell_idx,
            'neighbor_enrichment': enrichment,
        })

    results_df = pd.DataFrame(results_dict_list)
    results_df = results_df.groupby(['label']).mean().reset_index().drop(columns=['cell_index'])
    results_df.to_csv(output_file, sep="\t", index=False)


def evaluate_CITE_seq(input_adata_file, output_file, num_nn, num_pcs, protein_h5ad_file):
    print(f"input_adata_file: {input_adata_file}")
    print(f"output_file: {output_file}")
    print(f"num_nn: {num_nn}")
    print(f"num_pcs: {num_pcs}")
    print(f"protein_h5ad_file: {protein_h5ad_file}")
    # Rest of your function code goes here

    adata = sc.read_h5ad(input_adata_file)
    protein_adata = sc.read_h5ad(protein_h5ad_file)

    results_dict_list = []
    max_num_pcs = adata.obsm['X_pca'].shape[1]
    if num_pcs > max_num_pcs:
        raise ValueError("Not enough PCs to subset")
    X_pca = adata.obsm['X_pca'][:, :num_pcs]
    nbrs = NearestNeighbors(n_neighbors=num_nn, algorithm='brute', n_jobs=-1).fit(X_pca)
    _, knn_indices = nbrs.kneighbors(X_pca)


    results_dict_list = []
    for cell_idx in range(adata.shape[0]):
        neighbors = knn_indices[cell_idx]
        protein_distance_sum = 0
        for neighbor in neighbors:
            protein_distance_sum += np.linalg.norm(protein_adata.X[cell_idx] - protein_adata.X[neighbor]) # 2 norm by default
        result_dict = {
            'cell_index': cell_idx,
            'total_protein_neighbor_distance': protein_distance_sum,
        }
        results_dict_list.append(result_dict)

    results_df = pd.DataFrame(results_dict_list)
    results_df.to_csv(output_file, sep="\t", index=False)



from typing import Optional, Iterable, Tuple, Union
from scipy.sparse import issparse, csc_matrix, csr_matrix



def clr_normalize(adata: AnnData, inplace: bool = True, axis: int = 0) -> Union[None, AnnData]:
    """
    Taken from muon
    https://github.com/scverse/muon/blob/master/muon/_prot/preproc.py
    Apply the centered log ratio (CLR) transformation
    to normalize counts in adata.X.

    Args:
        data: AnnData object with protein expression counts.
        inplace: Whether to update adata.X inplace.
        axis: Axis across which CLR is performed.
    """

    if axis not in [0, 1]:
        raise ValueError("Invalid value for `axis` provided. Admissible options are `0` and `1`.")

    if not inplace:
        adata = adata.copy()

    if issparse(adata.X) and axis == 0 and not isinstance(adata.X, csc_matrix):
        x = csc_matrix(adata.X)
    elif issparse(adata.X) and axis == 1 and not isinstance(adata.X, csr_matrix):
        x = csr_matrix(adata.X)
    else:
        x = adata.X

    if issparse(x):
        x.data /= np.repeat(
            np.exp(np.log1p(x).sum(axis=axis).A / x.shape[axis]), x.getnnz(axis=axis)
        )
        np.log1p(x.data, out=x.data)
    else:
        np.log1p(
            x / np.exp(np.log1p(x).sum(axis=axis, keepdims=True) / x.shape[axis]),
            out=x,
        )

    adata.X = x

    return None if inplace else adata





