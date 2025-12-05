
import scanpy as sc
from sklearn import metrics
import numpy as np
import warnings
warnings.filterwarnings('ignore')
import squidpy as sq
#import SpaceFlow
import scanpy as sc
import pandas as pd

from SpaceFlow import SpaceFlow 
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
        for res in sorted(list(np.arange(0.2, 2.5, increment)), reverse=True):
            sc.tl.leiden(adata, random_state=0, resolution=res)
            count_unique_leiden = len(pd.DataFrame(adata.obs['leiden']).leiden.unique())
            if count_unique_leiden == fixed_clus_count:
                break
        return res



def run(adata ,data_name,data_type='Visium',n_clusters=7):
    start_time = time.time()
    tracemalloc.start()


    sf = SpaceFlow.SpaceFlow(adata=adata)

    #preprocess
    sf.preprocessing_data(n_top_genes=3000)

    """ sf.train(spatial_regularization_strength=0.1, 
         z_dim=50, 
         lr=1e-3, 
         epochs=1000, 
         max_patience=50, 
         min_stop=100, 
         random_seed=42, 
         #gpu=1, 
         regularization_acceleration=True, 
         edge_subset_sz=1000000) """
    adata.obsm['SpaceFlow_embedding']=sf.train()
    current, peak = tracemalloc.get_traced_memory()
    end_time = time.time()
    tracemalloc.stop()

    finaltime = f"{end_time - start_time:.4f}"
    current=f"{current / 10**6:.4f}"
    peak=f"{peak / 10**6:.4f}"
    print(f"Execution time: {finaltime} seconds")
    print(f"Current memory usage: {current} MB")
    print(f"Peak memory usage: {peak} MB")   
  
    sc.pp.neighbors(adata, n_neighbors=50)
    eval_resolution = res_search_fixed_clus(adata, n_clusters)

    sf.segmentation(domain_label_save_filepath="./dlpfc/dlpfc_domains_{}.csv", 
                n_neighbors=50, 
                resolution=eval_resolution)

    sf.plot_segmentation(segmentation_figure_save_filepath="./dlpfc/dlpfc_domain_segmentation_{}.pdf", 
                     colormap="tab20", 
                     scatter_sz=1., 
                     rsz=4., 
                     csz=4., 
                     wspace=.4, 
                     hspace=.5, 
                     left=0.125, 
                     right=0.9, 
                     bottom=0.1, 
                     top=0.9)
    
    pred=pd.read_csv('./dlpfc/dlpfc_domains_{}.csv',header=None)
    pred_list=pred.iloc[:,0].to_list()
    adata.obs['SpaceFlow_refine_domain'] = pred_list


    adata.obs['cluster'] =adata.obs['SpaceFlow_refine_domain'] 
    adata.uns['exec_time'] = finaltime
    adata.uns['current_memory'] = current   
    adata.uns['peak_memory'] = peak

    return adata.obs['SpaceFlow_refine_domain'],finaltime, peak
   
 