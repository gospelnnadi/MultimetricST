
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
Execute the following command to run all available spatial domain identification methods described in the paper or use --subset_methods to execute on a subset of methods.
`````
python MultimetricST.py \\<br>
  --mode 1 \\<br>
  --data_path Data/DLPFC/151673 \\<br>
  --ground_truth Data/DLPFC/151673/metadata.tsv \\<br>
  --ground_truth_col_name layer_guess \\<br>
  --data_name 151673 \\<br>
  --is_h5ad 0 \\<br>
  --data_type Visium \\<br>
  --n_clusters 7 &
`````

Example 2: Visium data stored as .h5ad. Download the Stereo Seq Axolotl Brain dataset found in the data availability section (see [README](README.md)).
Execute the following command to run all available spatial domain identification methods described in the paper or use `--subset_methods` to execute on a subset of methods. The parameter `--plot_size` is dataset specific and it required where the tissue image information is unavailable in adata.uns.  
`````
python MultimetricST.py \\<br>
  --mode 1 \\<br> 
  --data_path Data/Stereo/Stereo_Axolotl_Brain.h5ad \\<br>
  --ground_truth ground_truth \\<br>
  --is_h5ad 1 \\<br>
  --data_name Axolotl_Brain \\<br>
  --data_type Stereo \\<br>
  --plot_size 35  \\<br>
  --n_clusters 8 &
`````

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
`````
python MultimetricST.py \\<br>
  --mode 2 \\<br>
  --data_path Data/DLPFC/151673 \\<br>
  --is_h5ad 0 \\<br>
  --data_type Visium \\<br>
  --method_cluster_label CCST conST DeepST GIST GraphST SCAN-IT SEDR SpaceFlow SpaGCN STAGATE \\<br>
  --ground_truth Data/DLPFC/151673/metadata.tsv \\<br>
  --ground_truth_col_name layer_guess &
 `````
Example 4: Evaluate cluster labels from CSV file
`````
python MultimetricST.py \\<br>
  --mode 2 \\<br>
  --data_path Data/DLPFC/151673 \\<br>
  --is_h5ad 0 \\<br>
  --data_type Visium \\<br>
  --cluster_label_path multimetricST_outputs/clustering_labels.csv \\<br>
  --cluster_label_col_names CCST conST DeepST GIST GraphST SCAN-IT SEDR SpaceFlow SpaGCN STAGATE \\<br>
  --ground_truth Data/DLPFC/151673/metadata.tsv \\<br>
  --ground_truth_col_name layer_guess &
`````

### Mode 3 – Visualization only

Required

- CSV file with precomputed evaluation scores

Optional

- AnnData object or spatial matrix for spatial visualization


Mode 3 – Visualization Only Examples
Example 5: Visualize precomputed evaluation results

python MultimetricST.py \\<br>
  --mode 3 \\<br>
  --result_savepath multimetricST_outputs/clustering_results.csv  &

Example 6: Visualization with spatial plots
`````
python MultimetricST.py \\<br>
  --mode 3 \\<br>
  --data_path Data/DLPFC/151673 \\<br>
  --result_savepath multimetricST_outputs/clustering_results.csv \\<br>
  --cluster_label_path multimetricST_outputs/clustering_labels.csv \\<br>
  --cluster_label_col_names CCST conST DeepST GIST GraphST SCAN-IT SEDR SpaceFlow SpaGCN STAGATE  \\<br>
  --ground_truth Data/DLPFC/151673/metadata.tsv \\<br>
  --ground_truth_col_name layer_guess &
`````

# Notes

## Dashboard

The interactive dashboard:

- summarizes metrics by category,

- visualizes spatial clusters,

- highlights top-performing methods,

- compares computational efficiency (time & memory).

It is implemented using Panel and Plotly and is launched automatically at the end of the execution. The dashboard can be accessed via browser on http://localhost:5006/


#### Ground-Truth Annotations

Ground truth can be provided as:

- an adata.obs key, or

- an external CSV / TSV file

Example:

 - --ground_truth ground_truth, or
 
 - --ground_truth Data/DLPFC/151673/metadata.tsv
   --ground_truth_col_names layer_guess

If ground truth is not available, annotation-dependent metrics are skipped automatically.

## Output Files
All outputs are stored under:

  - multimetricST_outputs/clustering_results.csv

  - multimetricST_outputs/figures/<dataset_name>/<method_name>.png



## Dataset-Specific Parameters

The dataset in the data availability section (see [README](README.md)). Require number of cluster and plot size parameters: 

- **`n_cluster`** — Expected number of spatial domains (clusters). The number of cluster parameter is passed to each spatial domain identification method.
- **`plot_size`** — Point size used for spatial visualization in the dashboard (Use recommended plot sizes below or adjust when sizing leads to poor visibility)

---

### DLPFC (10x Visium Human Dorsolateral Prefrontal Cortex)

The number of expected clusters depends on the tissue section:

| Section IDs | `n_cluster` |`plot_size` |
|-------------|-------------|-------------|
| 151669, 151670, 151671, 151672 | 5 | 0 |
| 151507, 151508, 151509, 151510, 151673, 151674, 151675, 151676 | 7 | 0 |

---

### Other 10x Visium Human and Mouse Datasets

| Dataset | `n_cluster` | `plot_size` |
|---------|-------------|-------------|
| Human Breast Cancer | 20 | 0 |
| Human Ovarian Cancer | 8 | 0 |
| Human Lymph Node | 8 | 0 |
| Mouse Brain Anterior | 52 | 0 |
| Mouse Kidney | 7 | 0 |

---

### SlideSeq Mouse Datasets

| Dataset | `n_cluster` | `plot_size` |
|---------|-------------|-------------|
| Mouse Hippocampus | 14 | 35 |
| Mouse Olfactory Bulb | 7 | 35 |

---

### STARMaps Mouse Datasets

| Dataset | `n_cluster` | `plot_size` |
|---------|-------------|-------------|
| Mouse Visual Cortex | 7 | 250 |


---
### BaristaSeq Mouse Datasets

| Dataset | `n_cluster` | `plot_size` |
|---------|-------------|-------------|
| Mouse Primary Cortex | 6 | 15 |


---

### StereoSeq Axolotl and Mouse Datasets

| Dataset | `n_cluster` | `plot_size` |
|---------|-------------|-------------|
| Axolotl Brain | 16 | 35 |
| Mouse Brain (general) | 11 | 200 |

---

### StereoSeq MOSTA Embryo Datasets

MOSTA embryo sections require a very small plotting size due to their extremely high spatial resolution.

| Dataset | `n_cluster` | `plot_size` |
|---------|-------------|-------------|
| E9.5_E1S1 | 12 | 1 |
| E9.5_E2S1 | 14 | 1 |
| E10.5_E1S1 | 13 | 1 |
| E10.5_E2S1 | 18 | 1 |
| E11.5_E1S1 | 19 | 1 |
| E12.5_E1S1 | 23 | 1 |
| E12.5_E2S1 | 18 | 1 |
| E13.5_E1S1 | 19 | 1 |

---

### Xenium Datasets

| Dataset Type | `n_cluster` | `plot_size` |
|--------------|-------------|-------------|
| Mouse Brain Partial Coronal | 20 | 10 |

---

The `n_cluster` values are derived from the biological annotations when available in the datasets. Advanced users can set different parameters to explore alternative biological assumptions or benchmarking scenarios.
