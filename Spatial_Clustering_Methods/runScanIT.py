

import scanpy as sc
from sklearn import metrics

#pre-processing
import gc
import scanit
import torch
import random
import scanpy as sc
import pandas as pd
import anndata
import numpy as np
from scipy import sparse
from sklearn.metrics import normalized_mutual_info_score, adjusted_rand_score
from sklearn.cluster import SpectralClustering, KMeans
import matplotlib.pyplot as plt
#import stlearn as st
from pathlib import Path
 #SOMDE
from somde import SomNode
import time
import tracemalloc
# Start measuring time and memory
start_time = time.time()
tracemalloc.start()
def res_search_fixed_clus(adata, fixed_clus_count, increment=0.02):
    '''
        arg1(adata)[AnnData matrix]
        arg2(fixed_clus_count)[int]
        
        return:
            resolution[int]
    '''
    for res in sorted(list(np.arange(0.15, 5, increment)), reverse=True):
    # for res in sorted(list(np.arange(0.2, 2.5, increment)), reverse=True):
        sc.tl.leiden(adata, random_state=0, resolution=res)
        count_unique_leiden = len(pd.DataFrame(adata.obs['leiden']).leiden.unique())
        # print(res)
        # print(count_unique_leiden)
        if count_unique_leiden == fixed_clus_count:
            break
    return res

def run (adata ,data_name,data_type='Visium',n_clusters=7):
    start_time = time.time()
    tracemalloc.start()
   
    pts = adata.obsm['spatial']
    df_sp = pd.DataFrame(data=adata.X.toarray(), columns=list(adata.var_names))
    
    som = SomNode(pts, 5)

    ndf,ninfo = som.mtx(df_sp.T)
    nres = som.norm()
    result, SVnum =som.run()
    result.to_csv(f'./somde_result.csv')

    n_sv_genes = 3000
    sc.pp.normalize_total(adata)
    df_somde = pd.read_csv(f'./somde_result.csv')
    sv_genes = list( df_somde['g'].values[:n_sv_genes] )
    adata = adata[:, sv_genes]
    sc.pp.log1p(adata)
    sc.pp.scale(adata)

    scanit.tl.spatial_graph(adata, method='alpha shape', alpha_n_layer=2, knn_n_neighbors=5)
   
    scanit.tl.spatial_representation(adata, n_h=30, n_epoch=2000, lr=0.001, device='cpu', n_consensus=1, projection='mds', 
                python_seed=0, torch_seed=0, numpy_seed=0)
    
    current, peak = tracemalloc.get_traced_memory()
    end_time = time.time()
    tracemalloc.stop()
    finaltime = f"{end_time - start_time:.4f}"
    current=f"{current / 10**6:.4f}"
    peak=f"{peak / 10**6:.4f}"
    print(f"Execution time: {finaltime} seconds")
    print(f"Current memory usage: {current} MB")
    print(f"Peak memory usage: {peak} MB")

    sc.pp.neighbors(adata, use_rep='X_scanit', n_neighbors=15)
    eval_resolution = res_search_fixed_clus(adata, n_clusters)
    sc.tl.leiden(adata, key_added="scanit_leiden", resolution=eval_resolution)
   

    adata.obs['cluster'] =adata.obs['scanit_leiden'] 
    adata.uns['exec_time'] = finaltime
    adata.uns['current_memory'] = current   
    adata.uns['peak_memory'] = peak

 
    return adata.obs['scanit_leiden'],finaltime, peak
 