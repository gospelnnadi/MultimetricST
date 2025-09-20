
import subprocess
import os

# Read repo URLs from file (one per line)
with open("repos_git.txt", "r") as f:
    repos = [line.strip() for line in f if line.strip()]

for repo in repos:
    repo_name = os.path.splitext(os.path.basename(repo))[0]  # get folder name
    if os.path.exists(repo_name):
        print(f"Skipping {repo_name} (already exists)")
        continue

    print(f"Cloning {repo} ...")
    subprocess.run(["git", "clone", "--depth", "1", repo], check=True)



import sys, os

ROOT = os.path.dirname(os.getcwd()) #os.getcwd() #os.path.dirname(os.path.abspath(__file__))  # project root

# Add all method directories to sys.path
method_dirs = [
    
    "DeepST",
    #"GraphST",
    "GIST",
    "SEDR",
    "SpaGCN",
    "STAGATE_pyG",
    "SCAN-IT",
    "CCST",
    "conST",      # note: conST is lowercase
    "SpaceFlow"
]
current_dir = os.getcwd()
for d in method_dirs:
    path = os.path.join(current_dir, d)
    if path not in sys.path:
        sys.path.append(path)
sys.path.append(ROOT)

from runDeepST import run as runDeepST
from runGraphST import run as runGraphST
from runGIST import run as runGIST
from runSEDR import run as runSEDR
from runSpaGCN import run as runSpaGCN
from runSTAGATE import run as runSTAGATE 
from runScanIT import run as runScanIT
from runCCST import run as runCCST
from runConST import run as runConST
from runSpaceFlow import run as runSpaceFlow


from Silhouette_Spatial_Score.silhouette_spatial import silhouette_spatial_score
from scipy.sparse import issparse
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.decomposition import PCA
from multimetricST.evaluate import evaluate_cluster
from multimetricST.utils import *


def fromlayerstonumber (path):
  annotation_path= f"{path}/metadata.tsv"
  df_meta = pd.read_csv(annotation_path, sep='\t')
  df_meta_layer = df_meta['layer_guess']
  res=[]
  for sub in df_meta_layer.values:
    if sub == 'Layer1':
      res.append(str(sub).replace('Layer1', '1'))
    elif sub == 'Layer2':
      res.append(str(sub).replace('Layer2', '2'))
    elif sub == 'Layer3':
      res.append(str(sub).replace('Layer3', '3'))
    elif sub == 'Layer4':
      res.append(str(sub).replace('Layer4', '4'))
    elif sub == 'Layer5':
      res.append(str(sub).replace('Layer5', '5'))
    elif sub == 'Layer6':
      res.append(str(sub).replace('Layer6', '6'))
    elif sub == 'WM':
      res.append(str(sub).replace('WM', '7'))
    elif str(sub)=='nan' :
      res.append( res[-1]) ##nan
  return res


def evaluate(
    path, data_name, adata,
    is_h5ad=False,
    n_components=20,
    random_seed=35,
    is_visium=True,
):
    adata_raw = read_adata(f"{path}/{data_name}", is_h5ad=is_h5ad)
    adata_raw = preprocess(adata, n_components=n_components, random_seed=random_seed)

    pred = adata.obs['cluster']

    # Assign clusters (can be replaced with Leiden/Louvain/mclust etc.)
    ground_truth = fromlayerstonumber(f"{path}/{data_name}") 
    adata_raw.obs['ground_truth'] = ground_truth
    adata_raw=adata_raw[adata.obs_index,:]
    pca_matrix=adata_raw.obsm["X_pca"]



    scores=evaluate_cluster( adata_raw,pred, ground_truth, pca_matrix, is_visium=is_visium,verbose=True,decimal=4)
    return scores

def evaluate(
    path, data_name, adata,
    is_h5ad=False,
    n_components=20,
    random_seed=35,
    is_visium=True,
):
    # Read raw data
    adata_raw = read_adata(f"{path}/{data_name}", is_h5ad=is_h5ad)
    
    # Preprocess raw data (not adata!)
    adata_raw = preprocess(adata_raw, n_components=n_components, random_seed=random_seed)

    # Ground truth
    ground_truth = fromlayerstonumber(f"{path}/{data_name}")
    adata_raw.obs['ground_truth'] = ground_truth 
    
    # Keep only spots present in adata
    common_idx = adata_raw.obs_names.intersection(adata.obs_names)
    adata_raw = adata_raw[common_idx].copy()
    adata = adata[common_idx].copy()
    ground_truth =adata_raw.obs['ground_truth'] 
    # Clusters from adata
    pred = adata.obs['cluster']
    # PCA matrix
    pca_matrix = adata_raw.obsm["X_pca"]

    scores=evaluate_cluster( adata_raw,pred, ground_truth, pca_matrix, is_visium=is_visium,verbose=True,decimal=4)
    return scores



path=f"/home/accounts/personale/nndgpl46/MultimetricST/Data/DLPFC"
data_name="151673"

import pandas as pd

methods = [
    ("DeepST", runDeepST),
    ("GraphST", runGraphST),
    ("GIST", runGIST),
    ("SEDR", runSEDR),
    ("SpaGCN", runSpaGCN),
    ("STAGATE", runSTAGATE),
    ("ScanIT", runScanIT),
    ("CCST", runCCST),
    ("ConST", runConST),
    ("SpaceFlow", runSpaceFlow)
]

results = []

for method_name, method_func in methods:
    print(f"\n\nRunning {method_name}...\n\n\n\n\n")
    adata = method_func(path, data_name,n_clusters=7)
    
    scores = evaluate(path, data_name, adata)
    result = {
        "method": method_name,
        **scores,
        "exec_time": adata.uns.get('exec_time', None),
        "current_memory": adata.uns.get('current_memory', None),
        "peak_memory": adata.uns.get('peak_memory', None)
    }
    results.append(result)

results_df = pd.DataFrame(results)
results_df.to_csv("clustering_results1.csv", index=False)