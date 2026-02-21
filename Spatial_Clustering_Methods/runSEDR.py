import scanpy as sc
import pandas as pd
from sklearn import metrics
import torch

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')
import SEDR

from sklearn.decomposition import PCA


import time
import tracemalloc
# Start measuring time and memory
start_time = time.time()
tracemalloc.start()



def run(adata ,data_name,data_type='Visium',n_clusters=7):  
    start_time = time.time()
    tracemalloc.start()
    random_seed = 2023
    SEDR.fix_seed(random_seed)
    # gpu
    device = 'cuda:0' if torch.cuda.is_available() else 'cpu'

    adata.layers['count'] = adata.X.copy()
    sc.pp.filter_genes(adata, min_cells=50)
    sc.pp.filter_genes(adata, min_counts=10)
    sc.pp.normalize_total(adata, target_sum=1e6)
    sc.pp.highly_variable_genes(adata, flavor="seurat_v3", layer='count', n_top_genes=2000)
    adata = adata[:, adata.var['highly_variable'] == True]
    sc.pp.scale(adata)

    n_comp = min(200, adata.shape[1]-1)
    # sklearn PCA is used because PCA in scanpy is not stable.
    adata_X = PCA(n_components=n_comp, random_state=42).fit_transform(adata.X)
    adata.obsm['X_pca'] = adata_X
    graph_dict = SEDR.graph_construction(adata, 12)

    sedr_net = SEDR.Sedr(adata.obsm['X_pca'], graph_dict, mode='clustering', device=device)
    using_dec = True
    if using_dec:
     sedr_net.train_with_dec(N=1)
    else:
        sedr_net.train_without_dec(N=1)
    sedr_feat, _, _, _ = sedr_net.process()
    
    adata.obsm['SEDR'] = sedr_feat
    current, peak = tracemalloc.get_traced_memory()
    end_time = time.time()
    tracemalloc.stop()

    finaltime = f"{end_time - start_time:.4f}"
    current=f"{current / 10**6:.4f}"
    peak=f"{peak / 10**6:.4f}"
    print(f"Execution time: {finaltime} seconds")
    print(f"Current memory usage: {current} MB")
    print(f"Peak memory usage: {peak} MB")

    SEDR.mclust_R(adata, n_clusters, use_rep='SEDR', key_added='SEDR')

    adata.obs['cluster'] =adata.obs['SEDR'] 
    adata.uns['exec_time'] = finaltime
    adata.uns['current_memory'] = current   
    adata.uns['peak_memory'] = peak
    return adata.obs['SEDR'],finaltime, peak