
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

raw expression matrix

spatial coordinates

optional tissue image

optional ground-truth annotations



Mode 1 – Full Pipeline Examples
Example 1: DLPFC Visium dataset (as used in the paper)
python MultimetricST.py \
  --mode 1 \
  --data_path Data/DLPFC/151673 \
  --ground_truth Data/DLPFC/151673/metadata.tsv\
  --ground_truth_col_name layer_guess \
  --data_name 151673 \
  --n_clusters 7 &> outputs/multimetricst_mode1.log &


Example 2: Visium data stored as .h5ad
python MultimetricST.py \
  --mode 1 \
  --data_path Data/Stereo/Stereo_Axolotl_Brain.h5ad\
  --ground_truth ground_truth \
  --is_h5ad 1 \
  --data_name Axolotl_Brain \
  --n_clusters 8 > outputs/multimetricst_mode1.log &


This command:

clones all method repositories,

executes spatial domain identification methods,

evaluates clustering performance,

generates spatial plots,

launches the interactive dashboard.

###  Mode 2 – Evaluation + visualization

Required

AnnData object containing:

raw expression matrix

spatial coordinates

cluster labels

optional ground-truth annotations

OR

External cluster label files (CSV / TSV / NPY)


Mode 2 – Evaluation + Visualization Examples
Example 3: Evaluate cluster labels stored in adata.obs
python MultimetricST.py \
  --mode 2 \
  --data_path Data/DLPFC/151673 \
  --method_cluster_label GraphST SEDR SpaceFlow > outputs/multimetricst_mode2.log &
 
Example 4: Evaluate cluster labels from CSV file
python MultimetricST.py \
  --mode 2 \
  --data_path Data/DLPFC/151673 \
  --cluster_label_path clustering_labels.csv \
  --cluster_label_col_names GIST GraphST SEDR > outputs/multimetricst_mode2.log &


### Mode 3 – Visualization only

Required

CSV file with precomputed evaluation scores

Optional

AnnData object or spatial matrix for spatial visualization


Mode 3 – Visualization Only Examples
Example 5: Visualize precomputed evaluation results
python MultimetricST.py \
  --mode 3 \
  --result_savepath clustering_results.csv  > outputs/multimetricst_mode3.log &

Example 6: Visualization with spatial plots
python MultimetricST.py \
  --mode 3 \
  --data_path Data/DLPFC/151673 \
  --result_savepath clustering_results.csv \
  --cluster_label_path clustering_labels.csv \
  --cluster_label_col_names GraphST SEDR  > outputs/multimetricst_mode3.log &


## Dashboard

The interactive dashboard:

summarizes metrics by category,

visualizes spatial clusters,

highlights top-performing methods,

compares computational efficiency (time & memory).

It is implemented using Panel and Plotly and is launched automatically when applicable.


#### Ground-Truth Annotations

Ground truth can be provided as:

an adata.obs key, or

an external CSV / TSV file

Example:

--ground_truth ground_truth
--ground_truth_col_names layer_guess

If ground truth is not available, annotation-dependent metrics are skipped automatically.

## Output Files
All outputs are stored under:


multimetricST_outputs/clustering_results.csv
multimetricST_outputs/figures/<dataset_name>/<method_name>.png


### Notes on Method Execution

Method repositories are cloned automatically into Spatial_Clustering_Methods/

Runtime and memory usage are recorded for each method

New Python-based methods can be added by:

adding the GitHub URL to Spatial_Clustering_Methods/repos_git.txt

implementing a method-specific run function

### Known Compatibility Notes

#### SpaceFlow

The following change is required for compatibility:

File:

Spatial_Clustering_Methods/SpaceFlow/SpaceFlow/SpaceFlow.py line 132 

Replace:

 sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, flavor='cell_ranger', subset=True)


With:

sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, subset=True)


#### SEDR 
File:

Spatial_Clustering_Methods/SEDR/SEDR/clustering_func.py line 52 

Action:

to be commented in order to use the r-base in the conda environment

