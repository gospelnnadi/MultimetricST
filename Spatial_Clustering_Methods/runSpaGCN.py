import os,csv,re
import pandas as pd
import numpy as np
import scanpy as sc
import math
import SpaGCN as spg
from scipy.sparse import issparse
import random, torch
import warnings
warnings.filterwarnings("ignore")
import matplotlib.colors as clr
import matplotlib.pyplot as plt

#In order to read in image data, we need to install some package. Here we recommend package "opencv"
#inatll opencv in python
#!pip3 install opencv-python
#import cv2

import time
import tracemalloc
# Start measuring time and memory
start_time = time.time()
tracemalloc.start()


def run(adata ,data_name,data_type='Visium',n_clusters=7): 
    start_time = time.time()
    tracemalloc.start()
  
    """ spatial=pd.read_csv(f"{path}/{data_name}/spatial/tissue_positions_list.csv",sep=",",header=None,na_filter=False,index_col=0) 


    adata.obs["x1"]=spatial[1]
    adata.obs["x2"]=spatial[2]
    adata.obs["x3"]=spatial[3]
    adata.obs["x4"]=spatial[4]
    adata.obs["x5"]=spatial[5]
    #Select captured samples
    adata=adata[adata.obs["x1"]==1]
    adata.var_names=[i.upper() for i in list(adata.var_names)]
    adata.var["genename"]=adata.var.index.astype("str")

    #Set coordinates
    adata.obs["x_array"]=adata.obs["x2"]
    adata.obs["y_array"]=adata.obs["x3"]
    adata.obs["x_pixel"]=adata.obs["x4"]
    adata.obs["y_pixel"]=adata.obs["x5"]
    x_array=adata.obs["x_array"].tolist()
    y_array=adata.obs["y_array"].tolist()
    x_pixel=adata.obs["x_pixel"].tolist()
    y_pixel=adata.obs["y_pixel"].tolist()

    adata.var_names=adata.var.index.astype("str")
    adata.var_names_make_unique() """
    spg.prefilter_genes(adata,min_cells=3) # avoiding all genes are zeros
    spg.prefilter_specialgenes(adata)
    #Normalize and take log for UMI
    sc.pp.normalize_per_cell(adata)
    sc.pp.log1p(adata)

    adata.obs['x_pixel']=adata.obsm['spatial'][:,0]
    adata.obs['y_pixel']=adata.obsm['spatial'][:,1]

    x_pixel=adata.obs["x_pixel"].tolist()
    y_pixel=adata.obs["y_pixel"].tolist()

    #Calculate adjacent matrix
    s=1
    b=49

    #If histlogy image is not available, SpaGCN can calculate the adjacent matrix using the fnction below
    adj=spg.calculate_adj_matrix(x=x_pixel,y=y_pixel, histology=False)
    np.savetxt(f'./adj.csv', adj, delimiter=',')
    #spatial domain detection using SpaGCN

    #expression data preprocessing
    adj=np.loadtxt(f'./adj.csv', delimiter=',')

    #set hyper-parameters
    p=0.5 
    #Find the l value given p
    l=spg.search_l(p, adj, start=0.01, end=1000, tol=0.01, max_run=100)

    #Set seed
    r_seed=t_seed=n_seed=100
    #Seaech for suitable resolution
    res=spg.search_res(adata, adj, l, n_clusters, start=0.7, step=0.1, tol=5e-3, lr=0.05, max_epochs=20, r_seed=r_seed, t_seed=t_seed, n_seed=n_seed)

    #run SpaGCN
    clf=spg.SpaGCN()
    clf.set_l(l)
    #Set seed
    random.seed(r_seed)
    torch.manual_seed(t_seed)
    np.random.seed(n_seed)
    #Run
    clf.train(adata,adj,init_spa=True,init="louvain",res=res, tol=5e-3, lr=0.05, max_epochs=200)
    y_pred, prob=clf.predict()

    
    adata.obs["pred"]= y_pred
    adata.obs["pred"]=adata.obs["pred"].astype('category')

    current, peak = tracemalloc.get_traced_memory()
    end_time = time.time()
    tracemalloc.stop()

    finaltime = f"{end_time - start_time:.4f}"
    current=f"{current / 10**6:.4f}"
    peak=f"{peak / 10**6:.4f}"
    print(f"Execution time: {finaltime} seconds")
    print(f"Current memory usage: {current} MB")
    print(f"Peak memory usage: {peak} MB")

    adj_2d=spg.calculate_adj_matrix(x=x_pixel,y=y_pixel, histology=False)
    refined_pred=spg.refine(sample_id=adata.obs.index.tolist(), pred=adata.obs["pred"].tolist(), dis=adj_2d, shape="hexagon")
    adata.obs["refined_pred"]=refined_pred
    adata.obs["refined_pred"]=adata.obs["refined_pred"].astype('category')
    adata.obs['cluster'] =adata.obs['refined_pred'] 
    adata.uns['exec_time'] = finaltime
    adata.uns['current_memory'] = current   
    adata.uns['peak_memory'] = peak
    return adata.obs['refined_pred'],finaltime, peak


