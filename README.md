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

# setup

install requirements for python==3.10.0
````
pip install -r requirements.txt

````
install MultimetricST package
````
pip install git+https://github.com/InfOmics/MultimetricST.git
````


Create a [conda](https://www.anaconda.com/docs/getting-started/miniconda/install) environment
`````
conda create -n MMST2 python==3.10.0 r-base==4.3.1 seaborn=0.13.2  scipy==1.9.3 somoclu==1.7.5  gcc -y
conda activate MMST  

`````
install requirements 
````
pip install -r requirements.txt
````

install pytorch libraries
````
pip uninstall -y torch_sparse torch_scatter torch_cluster torch_spline_conv torch_geometric 

pip install torch==2.1.0 

pip install torch_sparse -f https://data.pyg.org/whl/torch-2.1.0+cu121.html
pip install torch_scatter -f https://data.pyg.org/whl/torch-2.1.0+cu121.html
pip install torch_cluster -f https://data.pyg.org/whl/torch-2.1.0+cu121.html
pip install torch_spline_conv -f https://data.pyg.org/whl/torch-2.1.0+cu121.html

pip install torch-geometric

````
Install mclust clustering algorithms used by some of the methods
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



### Data Availability ###
The spatial transcriptomics datasets are available at:  https://zenodo.org/records/17167458

Download the DLPFC 151673 data used in Clustering_Tutorial.py :
        wget https://zenodo.org/records/17167458/files/Data.zip

For detailed command-line usage and dataset-specific examples, see [View Usage Documentation](USAGE.md).