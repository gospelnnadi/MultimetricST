

import scanpy as sc
import numpy as np
import pandas as pd
from GIST.utils.clustering import cluster_n_plot
from GIST.GIST import GIST
import torch
import time
import tracemalloc



def read_adata(path, is_h5ad=False):
    #Todo: add check if Visium, if 'filtered_feature_bc_matrix.h5' or 'data_name_filtered_feature_bc_matrix.h5'
    if is_h5ad:
        adata = sc.read_h5ad(path)
        adata.var_names_make_unique()
        adata.obsm["spatial"]=adata.obsm["spatial"].astype(float)
    else: 
        adata = sc.read_visium(path, count_file='filtered_feature_bc_matrix.h5', load_images=True)
        adata.var_names_make_unique()
        adata.obsm["spatial"]=adata.obsm["spatial"].astype(float)
    return adata


def get_adata(path,  is_h5ad=False):


    if path=='':
        adata = read_adata('inputs/spatial_data/Data/1.DLPFC/151673' )
    else:
        adata =read_adata(path, is_h5ad)    
    return adata
        


def run(adata ,data_name,data_type='Visium',n_clusters=7):
 
    device =  "cuda" if torch.cuda.is_available() else "cpu"

    seed=35
    data_type=data_type
    refinement=True
    adata_raw=adata.copy()

    # Start measuring time and memory
    start_time = time.time()
    tracemalloc.start()
    

    GISTModel=GIST(adata=adata,  device=device, random_seed=seed, data_type=data_type)
    adata=GISTModel.train()

    current, peak = tracemalloc.get_traced_memory()
    end_time = time.time()
    tracemalloc.stop()

    finaltime = f"{end_time - start_time:.4f}"
    current=f"{current / 10**6:.4f}"
    peak=f"{peak / 10**6:.4f}"
    print(f"Execution time: {finaltime} seconds")
    print(f"Current memory usage: {current} MB")
    print(f"Peak memory usage: {peak} MB")  

    cluster_n_plot(adata, f"outputs/{data_name}.png", n_clusters,refinement=refinement, seed=seed, is_visium=GISTModel.is_visium)
    # Isolated spots are not considered in the clustering. If necessary copy the cluster label in the raw adata 
    adata_raw.obs['cluster']='-1'
    common = adata.obs_names.intersection(adata_raw.obs_names)
    adata_raw.obs.loc[common, 'cluster'] = adata.obs.loc[common, 'cluster'].values
    adata_raw.uns['GIST_emb']=adata.obsm['GIST_emb']
    
    adata.uns['exec_time'] = finaltime
    adata.uns['current_memory'] = current   
    adata.uns['peak_memory'] = peak
    return adata_raw.obs['cluster'],finaltime, peak
    

