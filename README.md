# MultimetricST Framework 
Multi-Metrics Framework for Spatial Transcriptomics clustering evaluation.

This repository contain metrics code from the article "Multi-Perspective Evaluation of Spatial Transcriptomics Clustering Methods"


install requirements for python==3.10.0
````
pip install -r requirements.txt

````
install MultimetricST package
````
pip install git+https://github.com/InfOmics/MultimetricST.git
````
test_mutlimetricST.py contains the MultimetricST Test pipeline on randomly generated cluster labels. 



### Given the predicted label and an optional.
    multimetricST/run_evaluator.py
    Evaluate clustering performance given an existing AnnData object and an optional raw dataset loading.
    Evaluate clustering performance given raw expression and spatial matrices.

### Spatial_Clustering_Methods folder contains the MultimetricST Framework for all 9 methods and the dashboard vissualizzation 
Create an environment if necessary
`````
conda create -n multimetricst python==3.10.0 r-base==4.3.1 -y

conda activate multimetricst

`````
install requirements 
````
pip install -r requirements2.txt

````

Run clustering and evaluation in folder Spatial_Clustering_Methods Clustering_Tutorial.py downloads the varous methods repos, run each method and  save results to clustering_results.csv:
    cd Spatial_Clustering_Methods

    python Clustering_Tutorial.py 
    
Visualize dashboard run ipynb script:
    runDashboard.ipynb

Note: To run the SpaceFlow method from the downloaded repo, the code In Spatial_Clustering_Methods/SpaceFlow/SpaceFlow/SpaceFlow.py line 132 need to use defualt flavour.  sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, flavor='cell_ranger', subset=True) -> sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, subset=True)

### Data Availability ###
The spatial transcriptomics datasets are available at:  https://zenodo.org/records/17167458

Download the DLPFC 151673 data used in Clustering_Tutorial.py :
        wget https://zenodo.org/records/17167458/files/Data.zip


