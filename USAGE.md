
---

# USAGE.md


# MultimetricST – Usage Guide

This document provides **practical usage instructions and example commands** for running MultimetricST in different modes and datasets, consistent with the experimental setup described in the paper.

---

## Execution Script

All functionality is accessed through:


python MultimetricST.py


## Input Data Requirements (by Mode)
###  Mode 1 – Full pipeline

Required

AnnData object containing:

- raw expression matrix

- spatial coordinates

- optional tissue image

- optional ground-truth annotations



Mode 1 – Full Pipeline Examples. 

Example 1: Download the DLPFC 10X Visium dataset found in the data availability section (see [README](README.md)).
Execute the following command.

python MultimetricST.py \\<br>
  --mode 1 \\<br>
  --data_path Data/DLPFC/151673 \\<br>
  --ground_truth Data/DLPFC/151673/metadata.tsv \\<br>
  --ground_truth_col_name layer_guess \\<br>
  --data_name 151673 \\<br>
  --is_h5ad 0 \\<br>
  --data_type Visium \\<br>
  --n_clusters 7 &


Example 2: Visium data stored as .h5ad. Download the Stereo Seq Axolotl Brain dataset found in the data availability section (see [README](README.md)).

python MultimetricST.py \\<br>
  --mode 1 \\<br> 
  --data_path Data/Stereo/Stereo_Axolotl_Brain.h5ad \\<br>
  --ground_truth ground_truth \\<br>
  --is_h5ad 1 \\<br>
  --data_name Axolotl_Brain \\<br>
  --data_type Stereo \\<br>
  --plot_size 35  \\<br>
  --n_clusters 8 &


This command:

clones all method repositories, executes spatial domain identification methods, evaluates clustering performance, generates spatial plots, launches the interactive dashboard.

###  Mode 2 – Evaluation + visualization

Required

AnnData object containing:

- raw expression matrix

- spatial coordinates

- cluster labels

- optional ground-truth annotations

OR

- External cluster label files (CSV / TSV / NPY)


Mode 2 – Evaluation + Visualization Examples.

Example 3: Evaluate cluster labels stored in adata.obs

python MultimetricST.py \\<br>
  --mode 2 \\<br>
  --data_path Data/DLPFC/151673 \\<br>
  --is_h5ad 0 \\<br>
  --data_type Visium \\<br>
  --method_cluster_label GraphST SEDR SpaceFlow &
 
Example 4: Evaluate cluster labels from CSV file

python MultimetricST.py \\<br>
  --mode 2 \\<br>
  --data_path Data/DLPFC/151673 \\<br>
  --is_h5ad 0 \\<br>
  --data_type Visium \\<br>
  --cluster_label_path clustering_labels.csv \\<br>
  --cluster_label_col_names GIST GraphST SEDR &


### Mode 3 – Visualization only

Required

- CSV file with precomputed evaluation scores

Optional

- AnnData object or spatial matrix for spatial visualization


Mode 3 – Visualization Only Examples
Example 5: Visualize precomputed evaluation results

python MultimetricST.py \\<br>
  --mode 3 \\<br>
  --result_savepath clustering_results.csv  &

Example 6: Visualization with spatial plots

python MultimetricST.py \\<br>
  --mode 3 \\<br>
  --data_path Data/DLPFC/151673 \\<br>
  --result_savepath clustering_results.csv \\<br>
  --cluster_label_path clustering_labels.csv \\<br>
  --cluster_label_col_names GraphST SEDR  \\<br>
  --ground_truth Data/DLPFC/151673/metadata.tsv \\<br>
  --ground_truth_col_name layer_guess &


## Dashboard

The interactive dashboard:

- summarizes metrics by category,

- visualizes spatial clusters,

- highlights top-performing methods,

- compares computational efficiency (time & memory).

It is implemented using Panel and Plotly and is launched automatically when applicable.


#### Ground-Truth Annotations

Ground truth can be provided as:

- an adata.obs key, or

- an external CSV / TSV file

Example:

--ground_truth ground_truth
--ground_truth_col_names layer_guess

If ground truth is not available, annotation-dependent metrics are skipped automatically.

## Output Files
All outputs are stored under:


multimetricST_outputs/clustering_results.csv
multimetricST_outputs/figures/<dataset_name>/<method_name>.png

