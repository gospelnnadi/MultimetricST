import time
import torch
# from memory_profiler import memory_usage
import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import sys
import torch.nn.functional as F
from sklearn.metrics.cluster import adjusted_rand_score
from sklearn.metrics import fowlkes_mallows_score
from sklearn.metrics import homogeneity_score
from sklearn.metrics import normalized_mutual_info_score
import math

import HERGAST

import time
import tracemalloc
# Start measuring time and memory
start_time = time.time()
tracemalloc.start()

def run(adata ,data_name,data_type='Visium',n_clusters=7): 
    start_time = time.time()
    tracemalloc.start()


    split = max(round(np.sqrt(adata.shape[0]/10000)),1)
    device_idx =  "cuda" if torch.cuda.is_available() else "cpu"

    ###preprocess
    n_comp = min(200, adata.shape[1]-1)
    sc.pp.normalize_total(adata, target_sum=1, exclude_highly_expressed=True)
    sc.pp.scale(adata)
    sc.pp.pca(adata, n_comps=n_comp)

    # save ground truth
    #sc.pl.embedding(adata, show=False, color='cell_type',basis='spatial',s=25,palette='tab20')
    #plt.savefig('GT_cell.pdf', bbox_inches='tight')

    # HERGAST pipeline
    HERGAST.utils.Cal_Spatial_Net(adata,verbose=False)
    HERGAST.utils.Cal_Expression_Net(adata, dim_reduce='PCA',verbose=False)

    train_HERGAST = HERGAST.Train_HERGAST(adata, batch_data=True, num_batch_x_y=(split,split), spatial_net_arg={'verbose':False},
                                      exp_net_arg={'verbose':False},dim_reduction='PCA',device_idx=device_idx)
    train_HERGAST.train_HERGAST(save_path=None, n_epochs=200,save_loss=False,save_reconstrction=False)

    current, peak = tracemalloc.get_traced_memory()
    end_time = time.time()
    tracemalloc.stop()

   
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

        # clustering
    res_list = [0.01,0.02,0.05,0.1,0.2,0.3,0.5,1]
    sc.pp.neighbors(adata, use_rep='HERGAST')
    res=0
    for used_res in res_list:
        res=used_res
        sc.tl.leiden(adata, random_state=2024, resolution=used_res,key_added='HERGAST')
        if n_clusters == len(adata.obs['HERGAST'].unique()):
            break

    print("Found resolution, ", res)
    adata.obs['cluster'] =adata.obs['HERGAST'] 
    adata.uns['exec_time'] = finaltime
    adata.uns['current_memory'] = current   
    adata.uns['peak_memory'] = peak
    return adata.obs['HERGAST'],finaltime, peak, allocated, cached


