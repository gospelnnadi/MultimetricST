


import os 
#from DeepST import run
import deepstkit as dt
import matplotlib.pyplot as plt
from pathlib import Path
import scanpy as sc


import time
import tracemalloc
# Start measuring time and memory
start_time = time.time()
tracemalloc.start()

def run(adata ,data_name,data_type='Visium',n_clusters=7,use_morphological=False):
    start_time = time.time()
    tracemalloc.start()
    SEED = 0   
    dt.utils_func.seed_torch(seed=SEED)

    #data_name = '2.Mouse_Brain_Anterior' #### project name
    save_path = "../Results" #### save path
    #n_clusters = 52 ###### the number of spatial domains.

    """  deepen = run(save_path = save_path,
	  task = "Identify_Domain", #### DeepST includes two tasks, one is "Identify_Domain" and the other is "Integration"
	  pre_epochs = 800, ####  choose the number of training
	  epochs = 1000, #### choose the number of training
	  use_gpu = False) """
    deepen = dt.main.run(
    #save_path=save_path,
    task="Identify_Domain",  # Spatial domain identification
    pre_epochs=500,          # Pretraining iterations
    epochs=500,              # Main training iterations
    use_gpu=True             # Accelerate with GPU if available
)
    ###### Read in 10x Visium data, or user can read in themselves.
    #adata = deepen._get_adata(platform="Visium", data_path=path, data_name=data_name)
    ###### Segment the Morphological Image
    #adata = deepen._get_image_crop(adata, data_name=data_name) 

    ###### Data augmentation. spatial_type includes three kinds of "KDTree", "BallTree" and "LinearRegress", among which "LinearRegress"
    ###### is only applicable to 10x visium and the remaining omics selects the other two.
    ###### "use_morphological" defines whether to use morphological images.
    #adata = deepen._get_augment(adata, spatial_type="LinearRegress", use_morphological=True)
    #adata = deepen._get_augment(adata, spatial_type="LinearRegress", use_morphological=False)

    
    if use_morphological:
    
      adata = deepen._get_image_crop(adata, data_name=data_name) 

      ###### Data augmentation. spatial_type includes three kinds of "KDTree", "BallTree" and "LinearRegress", among which "LinearRegress"
      ###### is only applicable to 10x visium and the remaining omics selects the other two.
      ###### "use_morphological" defines whether to use morphological images.
      adata = deepen._get_augment(adata, spatial_type="LinearRegress", use_morphological=True)
    else :
        adata.var_names_make_unique()
        ###### Data augmentation. spatial_type includes three kinds of "KDTree", "BallTree" and "LinearRegress", among which "LinearRegress"
        ###### is only applicable to 10x visium and the remaining omics selects the other two.
        ###### "use_morphological" defines whether to use morphological images.
        adata = deepen._get_augment(adata,use_data = "raw",  spatial_type = "BallTree", use_morphological=False)
    
    ###### Build graphs. "distType" includes "KDTree", "BallTree", "kneighbors_graph", "Radius", etc., see adj.py
    graph_dict = deepen._get_graph(adata.obsm["spatial"], distType = "BallTree")
  

    ###### Enhanced data preprocessing
    data = deepen._data_process(adata, pca_n_comps = 200)

    ###### Training models
    deepst_embed = deepen._fit(
		data = data,
		graph_dict = graph_dict,)
    ###### DeepST outputs
    adata.obsm["DeepST_embed"] = deepst_embed
    current, peak = tracemalloc.get_traced_memory()
    end_time = time.time()
    tracemalloc.stop()

    finaltime = f"{end_time - start_time:.4f}"
    current=f"{current / 10**6:.4f}"
    peak=f"{peak / 10**6:.4f}"
    print(f"Execution time: {finaltime} seconds")
    print(f"Current memory usage: {current} MB")
    print(f"Peak memory usage: {peak} MB")
    ###### Define the number of space domains, and the model can also be customized. If it is a model custom priori = False.
    adata = deepen._get_cluster_data(adata, n_domains=n_clusters, priori = True)

    ###### Spatial localization map of the spatial domain
    #sc.pl.spatial(adata, color='DeepST_refine_domain', frameon = False, spot_size=150)
    #plt.savefig(os.path.join(save_path, f'{data_name}_domains.pdf'), bbox_inches='tight', dpi=300)
    adata.obs['cluster'] =adata.obs['DeepST_refine_domain'] 
    adata.uns['exec_time'] = finaltime
    adata.uns['current_memory'] = current   
    adata.uns['peak_memory'] = peak
    return adata.obs['DeepST_refine_domain'],finaltime, peak

