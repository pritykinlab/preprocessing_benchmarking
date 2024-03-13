# import os
# import glob
# import pandas as pd
# import seaborn as sns
# import matplotlib.pyplot as plt

# # Base directory for the results and plots
# results_base_dir = "/Genomics/pritykinlab/yujie/preprocessing_benchmarking/results"
# plots_base_dir = os.path.join(results_base_dir, "plots")

# # Ensure the base plots directory exists
# os.makedirs(plots_base_dir, exist_ok=True)

# # Iterate over result files
# for result_file in glob.glob("/Genomics/pritykinlab/yujie/preprocessing_benchmarking/results/*/aggregated_results.tsv"):
#     dataset_name = result_file.split('/')[-2]  # Assumes a certain depth in the path
#     print(f"Processing: {dataset_name}")
    
#     # Read the results
#     results_df = pd.read_csv(result_file, sep='\t')

#     # Unique combinations for loop
#     labels = results_df['label'].unique()
#     hvg_norm_combos = results_df['hvg_norm_combo'].unique()
#     num_PCs_values = results_df['num_PCs'].unique()
#     num_nn_values = results_df['num_nn'].unique()

#     # Iterate over each combination and generate a plot
#     for label, hvg_norm_combo, num_PCs, num_nn in product(labels, hvg_norm_combos, num_PCs_values, num_nn_values):
#         # Filter the dataframe
#         filtered_df = results_df[(results_df['label'] == label) &
#                                  (results_df['hvg_norm_combo'] == hvg_norm_combo) &
#                                  (results_df['num_PCs'] == num_PCs) &
#                                  (results_df['num_nn'] == num_nn)]
        
#         if not filtered_df.empty:
#             # Plotting
#             plt.figure()
#             sns.lineplot(data=filtered_df, x='num_hvg', y='neighbors_count')
#             plt.title(f"{dataset_name} PCs: {num_PCs} K: {num_nn}")
#             plt.xlabel('Number of HVGs')
#             plt.ylabel('Neighbors Count')

#             # Directory for this specific plot
#             plot_dir = os.path.join(plots_base_dir, dataset_name, hvg_norm_combo, f"{num_PCs}_{num_nn}")
#             os.makedirs(plot_dir, exist_ok=True)

#             # Save the plot
#             plot_file_path = os.path.join(plot_dir, f"{label}.png")
#             plt.savefig(plot_file_path)
#             plt.close()


import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import glob

# Function to ensure directory exists
def ensure_dir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

# Base directory for the plots
results_base_dir = "/Genomics/pritykinlab/yujie/preprocessing_benchmarking/results"
base_plot_dir = os.path.join(results_base_dir, "plots")


# Iterate over each result file
for file_ in glob.glob("/Genomics/pritykinlab/yujie/preprocessing_benchmarking/results/*/aggregated_results.tsv"):
    dataset_name = file_.split('/')[-2]
    print(f"Processing: {dataset_name}")
    
    results_df = pd.read_csv(file_, sep='\t')

    # Iterate over unique combinations of hvg_norm_combo, num_PCs, num_nn
    for hvg_norm_combo in results_df['hvg_norm_combo'].unique():
        for num_PCs in results_df['num_PCs'].unique():
            for num_nn in results_df['num_nn'].unique():
                
                # Filter the dataframe for the current combination
                df_filtered = results_df[(results_df['hvg_norm_combo'] == hvg_norm_combo) &
                                         (results_df['num_PCs'] == num_PCs) &
                                         (results_df['num_nn'] == num_nn)]
                
                # Iterate over each label to generate a separate plot
                for label in df_filtered['label'].unique():
                    df_label = df_filtered[df_filtered['label'] == label]
                    
                    # Plotting
                    plt.figure()
                    sns.lineplot(x="num_hvg", y="neighbors_count", data=df_label)
                    plt.title(f'{dataset_name} PCs: {num_PCs} K: {num_nn}')
                    plt.xlabel('Number of HVGs')
                    plt.ylabel('Neighbors Count')

                    # Define the plot directory and filename
                    plot_dir = os.path.join(base_plot_dir, dataset_name, hvg_norm_combo, f"{num_PCs}_{num_nn}")
                    plot_filename = f"{label}.png"
                    ensure_dir(plot_dir)  # Ensure the directory exists

                    # Save the plot
                    plt.savefig(os.path.join(plot_dir, plot_filename))
                    plt.close()

