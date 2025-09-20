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

def run_evaluation(
    path='Data/DLPFC/151673',
    is_h5ad=False,
    n_components=20,
    random_seed=35,
    is_visium=True,
    n_clusters=7
):
    adata = read_adata(path, is_h5ad=is_h5ad)
    adata = preprocess(adata, n_components=n_components, random_seed=random_seed)

    # Assign clusters (can be replaced with Leiden/Louvain/mclust etc.)
    ground_truth = fromlayerstonumber(path) 
    pred = np.random.randint(0, n_clusters, size=adata.shape[0])
    pca_matrix=adata.obsm["X_pca"]

    evaluate_cluster( adata,pred, ground_truth, pca_matrix, is_visium=True,verbose=True,decimal=4)

    return 

run_evaluation()