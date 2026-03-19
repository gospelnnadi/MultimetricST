
import os
import torch
import ot
import pandas as pd
import anndata as ad
import numpy as np
import scanpy as sc
from sklearn import metrics
import multiprocessing as mp
from GraphST import GraphST
# clustering
from GraphST.utils import clustering


import time
import tracemalloc
# Start measuring time and memory


def run(adata ,data_name,data_type='Visium',n_clusters=7): 
    start_time = time.time()
    tracemalloc.start()
    device = 'cuda:0' if torch.cuda.is_available() else 'cpu'

    # define model
    model = GraphST.GraphST(adata, device=device,datatype=data_type)

    # train model
    adata = model.train()
    current, peak = tracemalloc.get_traced_memory()
    end_time = time.time()
    tracemalloc.stop()

    device_idx =  "cuda" if torch.cuda.is_available() else "cpu"
    if torch.cuda.is_available():
            device = torch.device(f"cuda:{device_idx}")
            allocated = torch.cuda.memory_allocated(device) / (1024 ** 2) 
            cached = torch.cuda.memory_reserved(device) / (1024 ** 2)
    else:
            allocated = cached = 0

    finaltime = f"{end_time - start_time:.4f}"
    current=f"{current / 10**6:.4f}"
    peak=f"{peak / 10**6:.4f}"
    print(f"Execution time: {finaltime} seconds")
    print(f"Current memory usage: {current} MB")
    print(f"Peak memory usage: {peak} MB")
    print(f"GPU memory allocated: {allocated:.4f} MB")
    print(f"GPU memory cached: {cached:.4f} MB")
    # set radius to specify the number of neighbors considered during refinement
    radius = 50

    tool = 'mclust' # mclust, leiden, and louvain
    print (f"Clustering using {tool}...")
    print("Dimension of Adata:",adata.obsm["emb"].shape)
    if tool == 'mclust':
        clustering(adata, n_clusters, radius=radius, method=tool, refinement=True) # For DLPFC dataset, we use optional refinement step.
    elif tool in ['leiden', 'louvain']:
        clustering(adata, n_clusters, radius=radius, method=tool, start=0.1, end=2.0, increment=0.01, refinement=False) 

    adata.obs['cluster'] =adata.obs['domain'] 
    adata.uns['exec_time'] = finaltime
    adata.uns['current_memory'] = current   
    adata.uns['peak_memory'] = peak
    return adata.obs['domain'],finaltime, peak, allocated, cached

