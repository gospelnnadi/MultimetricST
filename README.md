# MultimetricST Framework 
**MultimetricST** is a Python-based framework for the **evaluation of spatial domain identification methods** in spatial transcriptomics.

The framework enables **consistent, reproducible, and extensible benchmarking** by providing an integrated environment for:
- cloning or installing method repositories,
- executing spatial clustering methods,
- preprocessing spatial transcriptomics data,
- computing a comprehensive set of clustering evaluation metrics,
- visualizing results through an interactive dashboard.

MultimetricST accompanies the study:

> *A Multi-Perspective Evaluation Framework of Spatial Transcriptomics Clustering Methods*

---

## Overview

MultimetricST integrates all spatial domain identification methods evaluated in the study (see *Spatial domain identification methods* in the manuscript) into a unified framework.  
Users can easily **clone, execute, and compare methods** under a controlled environment, ensuring fair and reproducible evaluation.

The framework is organized into **three main components**:

1. **Method execution module**
2. **Cluster label evaluation module**
3. **Dashboard visualization module**

These components can be executed **independently or sequentially** via a single command-line:

# Setup
Setup notes: setup supported on Linux platform.

install MultimetricST package.
````
git clone https://github.com/InfOmics/MultimetricST.git
cd ~/MultimetricST
````

Download packages of the spatial transcriptomics spatial domain identification methods to be evaluated described in the paper.
````
python download_repo.py
````


Create a [conda](https://www.anaconda.com/docs/getting-started/miniconda/install#quickstart-install-instructions) (version>=23.11.0 is recommend) environment.
`````
conda create -n MMST -y

`````
`````
conda activate MMST  

`````
`````
conda install python=3.10.0 r-base=4.3.1 somoclu=1.7.5 -y 

`````

 ### install dependencies

````
pip install torch==2.1.0 

pip install torch_sparse -f https://data.pyg.org/whl/torch-2.1.0+cu121.html
pip install torch_scatter -f https://data.pyg.org/whl/torch-2.1.0+cu121.html
pip install torch_cluster -f https://data.pyg.org/whl/torch-2.1.0+cu121.html
pip install torch_spline_conv -f https://data.pyg.org/whl/torch-2.1.0+cu121.html

pip install torch-geometric==2.7.0 torchaudio==2.1.0 torchvision==0.16.0 
pip install -r requirements.txt
````

Install mclust clustering algorithms used by some of the spatial transcriptomics spatial domain identification methods.
````
Rscript -e 'install.packages("remotes", repos="https://cran.r-project.org")'

Rscript -e 'remotes::install_version("mclust", version = "6.0.1", repos="https://cran.r-project.org")'

````


## Execution Modes

The framework supports **three execution modes**, controlled by a user-specified parameter:

1. **Full pipeline** 
   Runs method execution (`clustering.py`), cluster evaluation (`evaluate.py`), and result visualization (`dashboard.py`).

2. **Evaluation + visualization**  
   Executes only cluster evaluation and dashboard visualization using precomputed cluster labels.

3. **Visualization only**  
   Visualizes results from a user-provided CSV file containing precomputed evaluation scores.

Each mode requires a different set of input data, consistent with the manuscript description.


### Notes on Method Execution

using the `download_repo.py` method repositories are cloned automatically into MultimetricST/Spatial_Clustering_Methods/

Runtime and memory usage are recorded for each method

New Python-based methods can be added by:

- adding the GitHub URL to Spatial_Clustering_Methods/repos_git.txt

- implementing a method-specific run function

### Known Compatibility Notes

The following change is required for compatibility. `download_repo.py` modifies the spatial domain identification method `SEDR` file: MultimetricST/Spatial_Clustering_Methods/SEDR/SEDR/clustering_func.py, line 52 is commented in order to use the r-base in the created conda environment (MMST). 



### Data Availability ###
The spatial transcriptomics datasets are available at:  https://zenodo.org/records/18482658

Download the DLPFC 151673 data:

        wget https://zenodo.org/records/18482658/files/Data.zip
        unzip Data.zip
        rm Data.zip

Download the Axolotl dataset:

        wget https://zenodo.org/records/18482658/files/Stereo.zip
        unzip Stereo.zip -d Data/
        rm Stereo.zip


For detailed command-line usage and dataset-specific examples, see [View Usage Documentation](USAGE.md).