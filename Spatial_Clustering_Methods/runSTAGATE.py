import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import os
import sys
import STAGATE_pyG 

import time
import tracemalloc
# Start measuring time and memory
start_time = time.time()
tracemalloc.start()

def run(adata ,data_name,data_type='Visium',n_clusters=7): 
    start_time = time.time()
    tracemalloc.start()


    #Normalization
    sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=3000)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    STAGATE_pyG .Cal_Spatial_Net(adata, rad_cutoff=150)
    STAGATE_pyG .Stats_Spatial_Net(adata)


    adata = STAGATE_pyG .train_STAGATE(adata)

    current, peak = tracemalloc.get_traced_memory()
    end_time = time.time()
    tracemalloc.stop()

    finaltime = f"{end_time - start_time:.4f}"
    current=f"{current / 10**6:.4f}"
    peak=f"{peak / 10**6:.4f}"
    print(f"Execution time: {finaltime} seconds")
    print(f"Current memory usage: {current} MB")
    print(f"Peak memory usage: {peak} MB")

    sc.pp.neighbors(adata, use_rep='STAGATE')
    sc.tl.umap(adata)
    adata = STAGATE_pyG.mclust_R(adata, used_obsm='STAGATE', num_cluster=n_clusters)
    adata.obs['cluster'] =adata.obs['mclust'] 
    adata.uns['exec_time'] = finaltime
    adata.uns['current_memory'] = current   
    adata.uns['peak_memory'] = peak
    return adata.obs['mclust'],finaltime, peak


