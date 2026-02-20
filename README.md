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

### Execution Modes

The framework supports **three execution modes**, controlled by a user-specified parameter:

1. **Full pipeline** 
   Runs method execution (`Spatial_Clustering_Methods/clustering.py`), cluster evaluation (`Evaluate/evaluate.py`), and result visualization (`Visualize_Scores/dashboard.py`).

2. **Evaluation + visualization**  
   Executes only cluster evaluation and dashboard visualization using precomputed cluster labels.

3. **Visualization only**  
   Visualizes results from a user-provided CSV file containing precomputed evaluation scores.
   

Each mode requires a different set of input data, consistent with the manuscript description.




# Setup

install MultimetricST package.
````
git clone https://github.com/InfOmics/MultimetricST.git
cd ~/MultimetricST
````

Download packages of the spatial transcriptomics spatial domain identification methods to be evaluated described in the paper.
````
python download_repo.py
````
If python not installed. 
````
sudo apt install python3-pip
````
````
python3 download_repo.py 
````
Repositories will be stored in:

MultimetricST/Spatial_Clustering_Methods/


## Setup - Conda 
The conda setup is supported on Linux platform. 

Create a [conda](https://www.anaconda.com/docs/getting-started/miniconda/install#quickstart-install-instructions) (version>=23.11.0 is recommend) environment. A conda quick install is also availabe in the install Miniconda section below.
`````
conda create -n MMST -y

`````
`````
conda activate MMST  

`````

 ### install dependencies

````
conda install -c conda-forge python=3.10.0 r-base=4.3.1 rpy2=3.5.11  r-essentials somoclu=1.7.5 -y 
pip install torch==2.1.0 

pip install torch_sparse -f https://data.pyg.org/whl/torch-2.1.0+cu121.html
pip install torch_scatter -f https://data.pyg.org/whl/torch-2.1.0+cu121.html
pip install torch_cluster -f https://data.pyg.org/whl/torch-2.1.0+cu121.html
pip install torch_spline_conv -f https://data.pyg.org/whl/torch-2.1.0+cu121.html

pip install torch-geometric==2.7.0 torchaudio==2.1.0 torchvision==0.16.0 
pip install -r requirements.txt

Rscript -e 'install.packages("remotes", repos="https://cran.r-project.org")'

Rscript -e 'remotes::install_version("mclust", version = "6.0.1", repos="https://cran.r-project.org")'

````


## Setup - Docker 
Setup on containerized environment using docker. To use the mmst container install [Docker](https://www.docker.com/) (version>=28.4.0 is recommend) and [Docker Compose] (https://docs.docker.com/compose/).

### Start the MultimetricST container
Get the mmst container with docker compose.
`````
docker-compose up -d
`````
Verify the container image.
````` 
docker images
`````


# Usage
##  Usage - Conda environment
Within the conda MMST environment (that is, conda activate MMST ). 

Download the DLPFC 10X Visium dataset found in the data availability section.
Execute the following command.

`````
python MultimetricST.py \
  --mode 1 \
  --data_path Data/DLPFC/151673 \
  --ground_truth Data/DLPFC/151673/metadata.tsv \
  --ground_truth_col_name layer_guess \
  --data_name 151673 \
  --is_h5ad 0 \
  --data_type Visium \
  --n_clusters 7 \
  --subset_methods CCST conST DeepST GIST GraphST SCAN-IT SEDR SpaceFlow SpaGCN STAGATE &
`````
At the end of the execution, the dashboard can be accessed via browser on http://localhost:5006/

For detailed command-line usage and dataset-specific and execution modes examples, see [View Usage Documentation](USAGE.md).

## Usage - Docker environment 
Access the container environment.
`````
docker-compose exec mmst bash
`````

Download the DLPFC 10X Visium dataset found in the data availability section.
Execute the following command.
`````
python MultimetricST.py \
  --mode 1 \
  --data_path Data/DLPFC/151673 \
  --ground_truth Data/DLPFC/151673/metadata.tsv \
  --ground_truth_col_name layer_guess \
  --data_name 151673 \
  --is_h5ad 0 \
  --data_type Visium \
  --n_clusters 7 \
  --subset_methods CCST conST DeepST GIST GraphST SCAN-IT SEDR SpaceFlow SpaGCN STAGATE &
  `````
At the end of the execution, the dashboard can be accessed via browser on http://localhost:5006/

For detailed command-line usage and dataset-specific and execution modes examples, see [View Usage Documentation](USAGE.md).

To exit the container environment. 
`````
Ctrl+D
`````
To stop the container.
`````
docker-compose down
`````



# Data Availability 
The spatial transcriptomics datasets are available at:  https://zenodo.org/records/18482658

Download the DLPFC 151673 data:

        wget https://zenodo.org/records/18482658/files/Data.zip
        unzip Data.zip
        rm Data.zip

Download the Axolotl dataset:

        wget https://zenodo.org/records/18482658/files/Stereo.zip
        unzip Stereo.zip -d Data/
        rm Stereo.zip


if unzip not installed.
````
sudo apt-get install zip unzip
````


# Notes
### Notes on Method Execution

using the `download_repo.py` method repositories are cloned automatically into MultimetricST/Spatial_Clustering_Methods/

Runtime and memory usage are recorded for each method

New Python-based methods can be added by:

- adding the GitHub URL to Spatial_Clustering_Methods/repos_git.txt

- Implementing a method-specific `run` function following the structure in  
  `Spatial_Clustering_Methods/run<MethodName>.py`

- Import the `run` function into `Spatial_Clustering_Methods/clustering.py` and add to the list of methods following the structure in `Spatial_Clustering_Methods/clustering.py`

All added spatial domain method repositories are executed through their respective `run` functions using **default parameters**.  

Users can customize and tune method-specific parameters by modifying the corresponding `run` function for each method.

`MultimetricST.py` allows users to specify dataset-specific parameters (e.g., **Visium**, **Slide**, **Stereo**, **Xenium**, **CosMx**, **VisiumHD**) through the `--data_type` argument. This parameter is passed to each method when required.

Using the `--subset_methods` parameter, users can define a subset of methods to be executed, evaluated, and visualized when running **Mode 1**.



### Known Compatibility Notes

The following modifications are automatically applied for compatibility within the MultimetricST environment:
   - **SEDR**  
         `download_repo.py` modifies the file `MultimetricST/Spatial_Clustering_Methods/SEDR/SEDR/clustering_func.py` by commenting out **line 52**. This change ensures compatibility with the **r-base** installation provided in the MMST Conda environment.

   - **SpaceFlow**  To resolve the error: "ValueError: Bin edges must be unique" which arises from the use of highly variable gene selection with `flavor='cell_ranger'`, the default parameter setting is used instead.  
 `download_repo.py` modifies the file `MultimetricST/Spatial_Clustering_Methods/SpaceFlow/SpaceFlow/SpaceFlow.py` (line 132).



### Information on the conda installation.
To ensure reproducibilty and easier setup we provide the conda quick installation

#### Install Miniconda
   ````
    mkdir -p ~/miniconda3
      wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
   bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
   rm ~/miniconda3/miniconda.sh
    ````

Make conda available
````
source ~/miniconda3/bin/activate
````
Verify conda installation
```` 
conda --version
```` 
#### Initialize conda
````
 conda init --all

 conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main
 conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r
 conda config --add channels conda-forge
````