
---
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




Mode 1 – Full Pipeline Examples.  at the end of the computation a **preprocessed AnnData object** is automatically saved in:
This allows **subsequent evaluation and export steps** without repeating the preprocessing phase, as the method cluster labels are stored directly in the saved AnnData object.

Example 1: Download the DLPFC 10X Visium dataset found in the data availability section (see [README](README.md)).
Execute the following command to run all available spatial domain identification methods described in the paper or use `--subset_methods` to execute on a subset of methods.
`````
python MultimetricST.py \
  --mode 1 \
  --data_path Data/DLPFC/151673 \
  --ground_truth Data/DLPFC/151673/metadata.tsv \
  --ground_truth_col_name layer_guess \
  --data_name DLPFC_151673 \
  --is_h5ad 0 \
  --result_filename clustering_results_DLPFC_151673.csv \
  --data_type Visium \
  --n_clusters 7 &
`````

Example 2: Download the Stereo Seq Axolotl Brain dataset found in the data availability section (see [README](README.md)).
Execute the following command to run all available spatial domain identification methods described in the paper or use `--subset_methods` to execute on a subset of methods. The parameter `--plot_size` is dataset specific and it required where the tissue image information is unavailable in adata.uns.  
`````
python MultimetricST.py \
  --mode 1 \ 
  --data_path Data/Stereo/Stereo_Axolotl_Brain.h5ad \
  --ground_truth ground_truth \
  --is_h5ad 1 \
  --data_name Axolotl_Brain \
  --result_filename clustering_results_Axolotl_Brain.csv \
  --data_type Stereo \
  --plot_size 35  \
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
python MultimetricST.py \
  --mode 2 \
  --data_path Data/Preprocessed/DLPFC_151673.h5ad \
  --is_h5ad 1 \
  --plot_size 0 \
  --data_name DLPFC_151673 \
  --data_type Visium \
  --result_filename clustering_results_DLPFC_151673.csv \
  --method_cluster_label CCST conST DeepST GIST GraphST HERGAST SCAN-IT SEDR SpaceFlow SpaGCN STAGATE \
  --ground_truth ground_truth  &
 `````
Example 4: Evaluate cluster labels from CSV file
`````
python MultimetricST.py \
  --mode 2 \
  --data_path Data/DLPFC/151673 \
  --is_h5ad 0 \
  --plot_size 0 \
  --data_name DLPFC_151673 \
  --data_type Visium \
  --result_filename clustering_results_DLPFC_151673.csv \
  --cluster_label_path multimetricST_outputs/clustering_labels.csv \
  --cluster_label_col_names CCST conST DeepST GIST GraphST HERGAST SCAN-IT SEDR SpaceFlow SpaGCN STAGATE \
  --ground_truth Data/DLPFC/151673/metadata.tsv \
  --ground_truth_col_name layer_guess &
`````

### Mode 3 – Visualization only

Required

- CSV file with precomputed evaluation scores

Optional

- AnnData object or spatial matrix for spatial visualization


Mode 3 – Visualization Only Examples

Example 5: Visualize precomputed evaluation results. 
`````
python MultimetricST.py \
  --mode 3 \
  --result_filename clustering_results_DLPFC_151673.csv \
  --data_name  DLPFC_151673 &
`````

Example 6: Visualization with spatial plots
`````
python MultimetricST.py \
  --mode 3 \
  --data_path Data/DLPFC/151673 \
   --is_h5ad 0 \
  --plot_size 0 \
  --result_filename clustering_results_DLPFC_151673.csv \
  --cluster_label_path multimetricST_outputs/clustering_labels.csv \
  --cluster_label_col_names CCST conST DeepST GIST GraphST HERGAST SCAN-IT SEDR SpaceFlow SpaGCN STAGATE  \
  --ground_truth Data/DLPFC/151673/metadata.tsv \
  --ground_truth_col_name layer_guess &
`````

# Notes


## Dashboard

The interactive dashboard:

- summarizes metrics by category,

- visualizes spatial clusters,

- highlights top-performing methods,

- compares computational efficiency (time & memory).

It is implemented using Panel and Plotly and is launched automatically at the end of the execution. The dashboard can be accessed via browser on http://localhost:8008/

Note: After successfully displaying the dashboard, launching a new dashboard instance may requires closing the browser tab of the previous instance (if any).

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
```
MultimetricST/multimetricST_outputs
```

Cluster evaluation score files are saved in:
  
  - multimetricST_outputs/clustering_results_<dataset_name>.csv

  filename specified using the `--result_filename` parameter

Cluster figures are saved in the directory:

  - multimetricST_outputs/figures/<dataset_name>/<method_name>.png

  filename is extracted automatically from `--data_name` parameter and methods.

For **mode 3**, the result file must be placed in the same directory (`MultimetricST/multimetricST_outputs`) to be correctly loaded.

If figures are available (or already saved) for the evaluated methods, for example in:

```
MultimetricST/multimetricST_outputs/figures/DLPFC_151673
```

they will be automatically displayed in the dashboard.

See `MultimetricST/run_simulated_datasets.ipynb` for an example workflow of **mode 3**.



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

### CosMx Datasets

| Dataset Type | `n_cluster` | `plot_size` |
|--------------|-------------|-------------|
| Human Lung Cancer | 18 | 100 |


---

The `n_cluster` values are derived from the biological annotations when available in the datasets. Advanced users can set different parameters to explore alternative biological assumptions or benchmarking scenarios.


## More Examples 
Mode 1 – Full Pipeline. 

Example 3: Download the DLPFC 10X Visium dataset found in the data availability section (see [README](README.md)).
Execute the following command to run all available spatial domain identification methods described in the paper or use `--subset_methods` to execute on a subset of methods.
`````
python MultimetricST.py \
  --mode 1 \
  --data_path Data/DLPFC/151674 \
  --ground_truth Data/DLPFC/151674/metadata.tsv \
  --ground_truth_col_name layer_guess \
  --data_name DLPFC_151674 \
  --is_h5ad 0 \
  --result_filename clustering_results_DLPFC_151674.csv \
  --data_type Visium \
  --n_clusters 7 &
`````

Example 4: Download the DLPFC 10X Visium dataset found in the data availability section (see [README](README.md)).
Execute the following command to run all available spatial domain identification methods described in the paper or use `--subset_methods` to execute on a subset of methods.
`````
python MultimetricST.py \
  --mode 1 \
  --data_path Data/DLPFC/151669 \
  --ground_truth Data/DLPFC/151669/metadata.tsv \
  --ground_truth_col_name layer_guess \
  --data_name DLPFC_151669 \
  --is_h5ad 0 \
  --result_filename clustering_results_DLPFC_151669.csv \
  --data_type Visium \
  --n_clusters 5 &
`````

Example 5: Download the DLPFC 10X Visium dataset found in the data availability section (see [README](README.md)).
Execute the following command to run all available spatial domain identification methods described in the paper or use `--subset_methods` to execute on a subset of methods.
`````
python MultimetricST.py \
  --mode 1 \
  --data_path Data/DLPFC/151508 \
  --ground_truth Data/DLPFC/151508/metadata.tsv \
  --ground_truth_col_name layer_guess \
  --data_name DLPFC_151508 \
  --is_h5ad 0 \
  --result_filename clustering_results_DLPFC_151508.csv \
  --data_type Visium \
  --n_clusters 7 &
`````



Example 6: Download the Mouse_Brain_Anterior 10X Visium dataset found in the data availability section (see [README](README.md)).
Execute the following command to run all available spatial domain identification methods described in the paper or use `--subset_methods` to execute on a subset of methods. 
`````
python MultimetricST.py \
  --mode 1 \ 
  --data_path Data/Mouse_Brain_Anterior \
  --ground_truth Data/Mouse_Brain_Anterior/metadata.tsv \
  --ground_truth_col_name ground_truth \
  --is_h5ad 0 \
  --data_name Mouse_Brain_Anterior \
  --data_type Visium \
  --result_filename clustering_results_Mouse_Brain_Anterior.csv \
  --plot_size 0  \
  --n_clusters 52 &
`````

Example 7: Download the STARmap dataset found in the data availability section (see [README](README.md)).
Execute the following command to run all available spatial domain identification methods described in the paper or use `--subset_methods` to execute on a subset of methods. The parameter `--plot_size` is dataset specific and it required where the tissue image information is unavailable in adata.uns.  
`````
python MultimetricST.py \
  --mode 1 \ 
  --data_path Data/STARmap_mouse_visual_cortex/STARmap_20180505_BY3_1k.h5ad \
  --ground_truth label \
  --is_h5ad 1 \
  --data_name Mouse_Visual_Cortex \
  --data_type STARmap \
  --result_filename clustering_results_Mouse_Visual_Cortex.csv \
  --plot_size 250  \
  --n_clusters 7 &
`````


Example 8: Download the BaristerSeq dataset found in the data availability section (see [README](README.md)).
Execute the following command to run all available spatial domain identification methods described in the paper or use `--subset_methods` to execute on a subset of methods. The parameter `--plot_size` is dataset specific and it required where the tissue image information is unavailable in adata.uns.  
`````
python MultimetricST.py \
  --mode 1 \ 
  --data_path Data/BaristerSeq-Mouse-Primary-Cortex/BaristerSeq_Mouse_Primary_Cortex.h5ad \
  --ground_truth ground_truth \
  --is_h5ad 1 \
  --data_name Mouse_Primary_Cortex \
  --data_type BaristerSeq \
  --result_filename clustering_results_BaristerSeq_Mouse_Primary_Cortex.csv \
  --plot_size 15  \
  --n_clusters 6 &
`````

Example 9: Download the MOSTA StereoSeq dataset found in the data availability section (see [README](README.md)).
Execute the following command to run all available spatial domain identification methods described in the paper or use `--subset_methods` to execute on a subset of methods. The parameter `--plot_size` is dataset specific and it required where the tissue image information is unavailable in adata.uns.  
`````
python MultimetricST.py \
  --mode 1 \ 
  --data_path Data/StereoSeq-MOSTA-Mouse_Organogenesis_Spatiotemporal_Transcriptomic_Atlas/E10.5_E2S1.MOSTA.h5ad \
  --ground_truth annotation \
  --is_h5ad 1 \
  --data_name MOSTA_E10_5_E2S1 \
  --data_type Stereo \
  --subset_methods conST DeepST GIST GraphST HERGAST SCAN-IT SEDR SpaceFlow SpaGCN \
  --result_filename clustering_results_MOSTA_E10_5_E2S1.csv \
  --plot_size 1  \
  --n_clusters 18 &
`````


Example 10: Download the MOSTA StereoSeq dataset found in the data availability section (see [README](README.md)).
Execute the following command to run all available spatial domain identification methods described in the paper or use `--subset_methods` to execute on a subset of methods. The parameter `--plot_size` is dataset specific and it required where the tissue image information is unavailable in adata.uns.  
`````
python MultimetricST.py \
  --mode 1 \ 
  --data_path Data/StereoSeq-MOSTA-Mouse_Organogenesis_Spatiotemporal_Transcriptomic_Atlas/E11.5_E1S1.MOSTA.h5ad \
  --ground_truth annotation \
  --is_h5ad 1 \
  --data_name MOSTA_E11_5_E1S1 \
  --subset_methods conST DeepST GIST GraphST HERGAST SCAN-IT SEDR SpaceFlow SpaGCN \
  --data_type Stereo \
  --result_filename clustering_results_MOSTA_E11_5_E1S1.csv \
  --plot_size 1  \
  --n_clusters 19 &
`````


Example 10: Download the Xenium Mouse Brain Partial Coronal dataset found in the data availability section (see [README](README.md)).
Execute the following command to run all available spatial domain identification methods described in the paper or use `--subset_methods` to execute on a subset of methods. The parameter `--plot_size` is dataset specific and it required where the tissue image information is unavailable in adata.uns.  
`````
python MultimetricST.py \
  --mode 1 \ 
  --data_path Data/Xenium/ms_brain_partialcoronal.h5ad \
  --is_h5ad 1 \
  --data_name Xenium_Mouse_Brain_Partial_Coronal \
  --data_type Xenium \
  --result_filename clustering_results_Xenium_Mouse_Brain_Partial_Coronal.csv \
  --plot_size 10  \
  --n_clusters 20 &
`````


Example 11: Download the CosMx Human Lung Cancer dataset found in the data availability section (see [README](README.md)).
Execute the following command to run all available spatial domain identification methods described in the paper or use `--subset_methods` to execute on a subset of methods. The parameter `--plot_size` is dataset specific and it required where the tissue image information is unavailable in adata.uns.  
`````
python MultimetricST.py \
  --mode 1 \ 
  --data_path Data/CosMx/Human_Lung_Cancer/SMI_Lung.h5ad \
   --ground_truth cell_type \
  --is_h5ad 1 \
  --data_name CosMx_Human_Lung_Cancer \
  --data_type CosMx \
  --subset_methods HERGAST SCAN-IT SEDR SpaceFlow STAGATE\
  --result_filename clustering_results_CosMx_Human_Lung_Cancer.csv \
  --plot_size 100  \
  --n_clusters 18 \
  --internal_metrics 0 &
`````


Example 12: Download the VisiumHD Human Colon Cancer dataset found in the data availability section (see [README](README.md)).
Execute the following command to run all available spatial domain identification methods described in the paper or use `--subset_methods` to execute on a subset of methods. The parameter `--plot_size` is dataset specific and it required where the tissue image information is unavailable in adata.uns.  
`````
python MultimetricST.py \
  --mode 1 \ 
  --data_path Data/HDVisium/Human_Colorectal_Cancer/binned_outputs/square_016um  \
  --is_h5ad 0 \
  --ground_truth Data/HDVisium/Human_Colorectal_Cancer/binned_outputs/16um_squares_annotation.csv \
  --data_name VisiumHD016_Human_Lung_Cancer \
  --data_type VisiumHD \
  --result_filename clustering_results_VisiumHD016_Human_Colon_Cancer.csv \
  --plot_size 100  \
  --n_clusters 18 \
   --internal_metrics 0 &
`````



## More Examples 
Mode 2 – Evaluation + visualization. --data_path if provided by the user it is the path to the adata containing the cluster labels and spatial location for optional spatial plotting visualizzation. The --method_cluster_label  are the name of the adata.obs key containing the cluster the key should correspond to the method name. a Prepocessed adata data in Data/Preprocessed is saved when the mode 1 is previuosly used on the same data set this allows subsequent evalution and export with method cluster labels saved. In the case a csv file contiain the clusters is available instead of the adata use the --cluster_label_path --cluster_label_col_names  option. 

## More Examples - Mode 2 – Evaluation + Visualization

In **Mode 2**, the framework performs **evaluation and optional visualization** of clustering results.
- `--data_path`:  
  If provided by the user, this should be the path to the **AnnData (`.h5ad`) object** containing the **cluster labels** and **spatial coordinates**, which are required for optional spatial visualization.
- `--method_cluster_label`:  
  These are the **names of the keys in `adata.obs` containing the cluster labels**.  
  Each key should correspond to the **name of the method** used to generate the clustering.
When **Mode 1** has previously been executed on the same dataset, a **preprocessed AnnData object** is automatically saved in:
This allows **subsequent evaluation and export steps** without repeating the preprocessing phase, as the method cluster labels can be stored directly in the saved AnnData object.
Alternatively, if the cluster assignments are stored in a **CSV file instead of an AnnData object**, you can use the following options:
- `--cluster_label_path` – path to the CSV file containing cluster labels  
- `--cluster_label_col_names` – names of the columns containing the cluster labels for each method



`````
python MultimetricST.py \
  --mode 2 \
  --data_path Data/Preprocessed/DLPFC_151673.h5ad \
  --is_h5ad 1 \
  --plot_size 0 \
  --data_name DLPFC_151673 \
  --data_type Visium \
  --result_filename clustering_results_evaluate_DLPFC_151673.csv \
  --method_cluster_label CCST conST DeepST GIST GraphST HERGAST SCAN-IT SEDR SpaceFlow SpaGCN STAGATE \
  --ground_truth ground_truth &
`````


`````
python MultimetricST.py \
  --mode 2 \
  --data_path Data/Preprocessed/DLPFC_151674.h5ad \
  --is_h5ad 1 \
  --plot_size 0 \
  --data_name DLPFC_151674 \
  --data_type Visium \
  --result_filename clustering_results_evaluate_DLPFC_151674.csv \
  --method_cluster_label CCST conST DeepST GIST GraphST HERGAST SCAN-IT SEDR SpaceFlow SpaGCN STAGATE \
  --ground_truth ground_truth &

`````


`````
python MultimetricST.py \
  --mode 2 \
  --data_path Data/Preprocessed/DLPFC_151508.h5ad \
  --is_h5ad 1 \
  --plot_size 0 \
  --data_name DLPFC_151508 \
  --data_type Visium \
  --result_filename clustering_results_evaluate_DLPFC_151508.csv \
  --method_cluster_label CCST conST DeepST GIST GraphST HERGAST SCAN-IT SEDR SpaceFlow SpaGCN STAGATE \
  --ground_truth ground_truth  &

`````


`````
python MultimetricST.py \
  --mode 2 \
  --data_path Data/Preprocessed/DLPFC_151669.h5ad \
  --is_h5ad 1 \
  --plot_size 0 \
  --data_name DLPFC_151669 \
  --data_type Visium \
  --result_filename clustering_results_evaluate_DLPFC_151669.csv \
  --method_cluster_label CCST conST DeepST GIST GraphST HERGAST SCAN-IT SEDR SpaceFlow SpaGCN STAGATE \
  --ground_truth ground_truth   &
`````


`````
python MultimetricST.py \
  --mode 2 \
  --data_path Data/Preprocessed/Mouse_Brain_Anterior.h5ad \
  --is_h5ad 1 \
  --plot_size 0 \
  --data_name Mouse_Brain_Anterior \
  --data_type Visium \
  --result_filename clustering_results_evaluate_Mouse_Brain_Anterior.csv \
  --method_cluster_label CCST conST DeepST GIST GraphST HERGAST SCAN-IT SEDR SpaceFlow SpaGCN STAGATE \
  --ground_truth ground_truth  &
`````

`````
python MultimetricST.py \
  --mode 2 \
  --data_path Data/Preprocessed/Mouse_Primary_Cortex.h5ad \
  --is_h5ad 1 \
  --plot_size 15 \
  --data_name Mouse_Primary_Cortex \
  --data_type BaristerSeq \
  --result_filename clustering_results_evaluate_Mouse_Primary_Cortex.csv \
  --method_cluster_label CCST conST DeepST GIST GraphST HERGAST SCAN-IT SEDR SpaceFlow SpaGCN STAGATE \
  --ground_truth ground_truth  &
`````

`````
python MultimetricST.py \
  --mode 2 \
  --data_path Data/Preprocessed/Mouse_Visual_Cortex.h5ad \
  --is_h5ad 1 \
  --plot_size 250 \
  --data_name Mouse_Visual_Cortex \
  --data_type STARMaps \
  --result_filename clustering_results_evaluate_Mouse_Visual_Cortex.csv \
  --method_cluster_label CCST conST DeepST GIST GraphST HERGAST SCAN-IT SEDR SpaceFlow SpaGCN STAGATE \
  --ground_truth ground_truth   &
`````




`````
python MultimetricST.py \
  --mode 2 \
  --data_path Data/Preprocessed/Axolotl_Brain.h5ad \
  --is_h5ad 1 \
  --plot_size 35 \
  --data_name Axolotl_Brain \
  --data_type Stereo \
  --result_filename clustering_results_evaluate_Axolotl_Brain.csv \
  --method_cluster_label CCST conST DeepST GIST GraphST HERGAST SCAN-IT SEDR SpaceFlow SpaGCN STAGATE \
  --ground_truth ground_truth  &
`````

`````
python MultimetricST.py \
  --mode 2 \
  --data_path Data/Preprocessed/MOSTA_E10_5_E2S1.h5ad\
  --is_h5ad 1 \
  --plot_size 1 \
  --data_name MOSTA_E10_5_E2S1 \
  --data_type Stereo \
  --result_filename clustering_results_evaluate_MOSTA_E10_5_E2S1.csv \
  --method_cluster_label conST DeepST GIST GraphST HERGAST SCAN-IT SEDR SpaceFlow SpaGCN \
  --ground_truth ground_truth   &
`````


`````
python MultimetricST.py \
  --mode 2 \
  --data_path Data/Preprocessed/Xenium_Mouse_Brain_Partial_Coronal.h5ad\
  --is_h5ad 1 \
  --plot_size 10 \
  --data_name Xenium_Mouse_Partial_Coronal \
  --data_type Xenium \
  --result_filename clustering_results_evaluate_Xenium_Mouse_Partial_Coronal.csv \
  --method_cluster_label CCST conST DeepST GIST GraphST HERGAST SCAN-IT SEDR SpaceFlow SpaGCN STAGATE \
  --ground_truth ground_truth   &

`````


`````
python MultimetricST.py \
  --mode 2 \
  --data_path Data/Preprocessed/MOSTA_E11_5_E1S1.h5ad\
  --is_h5ad 1 \
  --plot_size 1 \
  --data_name MOSTA_E11_5_E1S1 \
  --data_type Stereo \
  --result_filename clustering_results_evaluate_MOSTA_E11_5_E1S1.csv \
  --method_cluster_label conST DeepST GIST GraphST HERGAST SCAN-IT SEDR SpaceFlow SpaGCN \
  --ground_truth ground_truth  &

`````


`````
python MultimetricST.py \
  --mode 2 \
  --data_path Data/Preprocessed/CosMx_Human_Lung_Cancerl.h5ad\
  --is_h5ad 1 \
  --plot_size 100 \
  --data_name CosMx_Human_Lung_Cancer \
  --data_type CosMx \
  --result_filename clustering_results_evaluate_CosMx_Human_Lung_Cancer.csv \
  --method_cluster_label HERGAST SCAN-IT SEDR SpaceFlow STAGATE \
  --ground_truth ground_truth\
  --internal_metrics 0  &

`````

`````
python MultimetricST.py \
  --mode 2 \
  --data_path Data/Preprocessed/VisiumHD016_Human_Colon_Cancer.h5ad\
  --is_h5ad 1 \
  --plot_size 0 \
  --data_name VisiumHD016_Human_Colon_Cancer \
  --data_type VisiumHD \
  --result_filename clustering_results_evaluate_VisiumHD016_Human_Colon_Cancer.csv \
  --method_cluster_label HERGAST SpaceFlow \
  --ground_truth ground_truth \
  --internal_metrics 0   &
`````