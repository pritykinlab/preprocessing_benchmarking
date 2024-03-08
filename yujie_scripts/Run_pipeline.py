import os
import scanpy as sc
import pandas as pd
import numpy as np
import scipy
from scipy.sparse import csr_matrix
from itertools import product
from sklearn.neighbors import NearestNeighbors

def run_pipeline(input_adata_file, output_dir, pipeline_params):
    """
    Executes a flexible pipeline on scRNA-seq data based on provided steps and parameters.
    """
    # Initial setup - preparing 'all_labels.txt' and 'cleaned_input.h5ad'
    prepare_all_labels(input_adata_file, output_dir, pipeline_params[-1][1]['label_col'][0])
    cleaned_input_adata_file = prepare_cleaned_input(input_adata_file, output_dir)
    
    all_output_adata_files_in_prev_step = [cleaned_input_adata_file]
    for step_index, (func, params) in enumerate(pipeline_params[:-1]):  # Exclude last step for separate processing
    # for func, params in pipeline_params[:-1]:  # the index (step_index in this case) is not needed for any specific operation within the loop, you can safely remove it along with the enumerate function.
        all_output_adata_files_in_cur_step = []
        for input_adata_file in all_output_adata_files_in_prev_step:
            base_file_name = os.path.basename(input_adata_file).replace('.h5ad', '')
            for comb_params in generate_combinations(params):
                output_file_name = f"{base_file_name}__" + "__".join([f"{key}%{value}" for key, value in comb_params.items()]) + ".h5ad"

                output_adata_file = os.path.join(output_dir, output_file_name)
                func(input_adata_file, output_adata_file, **comb_params)
                all_output_adata_files_in_cur_step.append(output_adata_file)
        all_output_adata_files_in_prev_step = all_output_adata_files_in_cur_step

    # Process the last step separately and save individual .tsv files
    all_metrics_dfs = []
    for adata_file in all_output_adata_files_in_prev_step:
        print(adata_file)
        func, params = pipeline_params[-1]
        for comb_params in generate_combinations(params):
            print(comb_params)
            # Construct a filename for the individual results based on parameters
            individual_tsv_filename = construct_individual_filename(adata_file, comb_params) + ".tsv"
            individual_tsv_path = os.path.join(output_dir, individual_tsv_filename)

            # Check if the individual .tsv file already exists
            if not os.path.exists(individual_tsv_path):
                metric_df = func(adata_file, **comb_params)
                for key, value in comb_params.items():
                    metric_df[key] = value
                params_from_file = extract_params_from_filename(adata_file)
                for key, value in params_from_file.items():
                    metric_df[key] = value
                # Save the individual DataFrame to a .tsv file
                metric_df.to_csv(individual_tsv_path, sep="\t", index=False)
            else:
                print(f"{individual_tsv_path} already exists. Skipping generation.")
                # Read the existing .tsv file into a DataFrame
                metric_df = pd.read_csv(individual_tsv_path, sep="\t")

            all_metrics_dfs.append(metric_df)

    # Aggregate all individual .tsv files into "aggregated_results.tsv"
    aggregated_df = pd.concat(all_metrics_dfs, ignore_index=True)
    aggregated_filename = os.path.join(output_dir, "aggregated_results.tsv")
    aggregated_df.to_csv(aggregated_filename, sep="\t", index=False)
    return aggregated_df

def construct_individual_filename(adata_file, comb_params):
    """
    Construct a unique filename for the individual .tsv files based on the input parameters.
    """
    base_name = os.path.basename(adata_file).replace('.h5ad', '')
    param_str = "__".join([f"{key}%{value}" for key, value in comb_params.items()])
    return f"{base_name}__{param_str}"


def generate_combinations(params):
    """Generate all combinations of parameters."""
    keys, values = zip(*params.items())
    return [dict(zip(keys, v)) for v in product(*values)]

def extract_params_from_filename(filename):
    """Extract parameters and their values from a filename."""
    # Initialize an empty dictionary to store extracted parameters
    params = {}

    # Remove the file extension and any directory paths
    base_filename = os.path.splitext(os.path.basename(filename))[0]

    # Split the base filename into key-value pairs on the double underscore '__'
    kv_pairs = base_filename.split('__')

    # For each key-value pair, split on the percent '%' to separate the key from the value
    for kv in kv_pairs:
        if '%' in kv:
            key, value = kv.split('%', 1)  # Only split on the first percent
            params[key] = value

    return params

def prepare_all_labels(input_adata_file, output_dir, label_col):
    """Prepare 'all_labels.txt' containing unique labels from the .obs attribute of the AnnData object."""
    
    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)  # This line ensures the directory is created if it does not exist
    
    all_labels_file = os.path.join(output_dir, 'all_labels.txt')
    if not os.path.exists(all_labels_file):
        adata = sc.read_h5ad(input_adata_file)      
        # Ensure label_col is a single column and extract unique values correctly
        all_labels = adata.obs[label_col].unique()
        
        with open(all_labels_file, 'w') as f:
            for label in all_labels:
                f.write(f"{label}\n")
    print("Finished preparing 'all_labels.txt'")

def prepare_cleaned_input(input_adata_file, output_dir):
    """Prepare 'cleaned_input.h5ad' from the input AnnData object."""
    output_adata_file = os.path.join(output_dir, "cleaned_input.h5ad")
    if not os.path.exists(output_adata_file):
        adata = sc.read_h5ad(input_adata_file)
        if not isinstance(adata.X, scipy.sparse.spmatrix):
            adata.X = csr_matrix(adata.X)
        adata.layers['raw'] = adata.X
        adata.write_h5ad(output_adata_file)
        print("Finished preparing 'cleaned_input.h5ad'")
    return output_adata_file


############ Function implementations for pipeline steps. Add your customized functions below ############

def hvg(input_adata_file, output_adata_file, hvg_method, num_hvg):
    # Check if the output file already exists
    if os.path.exists(output_adata_file):
        print(f"File {output_adata_file} already exists. Skipping HVG.")
        return
        
    adata = sc.read_h5ad(input_adata_file)
    sc.pp.filter_genes(adata, min_cells=3)

    adata_original = adata.copy()
    # Apply HVG selection based on the method
    if hvg_method in ['seurat', 'cell_ranger', 'seurat_v3', 'pearson_residuals', 'pearson_residuals_1']:
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
        elif hvg_method == 'pearson_residuals' or hvg_method == 'pearson_residuals_1':
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



def evaluate(input_adata_file, label_col, num_nn, num_pcs_list):

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

    print("Finished evaluation")
    return results_df


# Example pipeline configuration
test_params= [
    # Specify each step with its parameters
    (hvg, {'hvg_method': ['seurat', 'cell_ranger'], 'num_hvg': [2000, 3000]}),
    (norm, {'norm_method': ['log_zscore', 'pearson_residuals']}),
    (pca, {'max_pcs': [50]}),
    (evaluate, {'label_col': ['gene'],'num_nn': [20],'num_pcs_list': [10, 25]})
]

# Uncomment the following line to run the pipeline with your specific input file and desired output directory
run_pipeline(input_adata_file="/Genomics/pritykinlab/share/perturbseq/GW_perturbseq/K562_essential_my_raw_all_genes_small_v2.h5ad", output_dir="/Genomics/pritykinlab/yujie/preprocessing_benchmarking/results4", pipeline_params=test_params)
