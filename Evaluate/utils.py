
import numpy as np
from scipy.spatial import *
from sklearn.preprocessing import *
from sklearn.metrics import *
from scipy.spatial.distance import *
from scipy.spatial import distance_matrix
import scipy.sparse as sp
from scipy.spatial.distance import cdist


from scipy.sparse import issparse
import scanpy as sc
from sklearn.decomposition import PCA


def compute_CHAOS(clusterlabel, location):

        clusterlabel = np.array(clusterlabel)
        location = np.array(location)
        matched_location = StandardScaler().fit_transform(location)

        clusterlabel_unique = np.unique(clusterlabel)
        dist_val = np.zeros(len(clusterlabel_unique))
        count = 0
        for k in clusterlabel_unique:
            location_cluster = matched_location[clusterlabel==k,:]
            if len(location_cluster)<=2:
                continue
            n_location_cluster = len(location_cluster)
            results = [fx_1NN(i,location_cluster) for i in range(n_location_cluster)]
            dist_val[count] = np.sum(results)
            count = count + 1

        return np.sum(dist_val)/len(clusterlabel)
    


def fx_1NN(i,location_in):
        location_in = np.array(location_in)
        dist_array = distance_matrix(location_in[i,:][None,:],location_in)[0,:]
        dist_array[i] = np.inf
        return np.min(dist_array)
    

def fx_kNN(i,location_in,k,cluster_in):

        location_in = np.array(location_in)
        cluster_in = np.array(cluster_in)


        dist_array = distance_matrix(location_in[i,:][None,:],location_in)[0,:]
        dist_array[i] = np.inf
        ind = np.argsort(dist_array)[:k]
        cluster_use = np.array(cluster_in)
        if np.sum(cluster_use[ind]!=cluster_in[i])>(k/2):
            return 1
        else:
            return 0
        

           
def compute_PAS(clusterlabel,location):
        
        clusterlabel = np.array(clusterlabel)
        location = np.array(location)
        matched_location = location
        results = [fx_kNN(i,matched_location,k=10,cluster_in=clusterlabel) for i in range(matched_location.shape[0])]
        return np.sum(results)/len(clusterlabel)
    
def compute_ASW(clusterlabel,location):
        return silhouette_score(X=location,labels=clusterlabel,metric='euclidean')  



def pca (X, n_components=20,random_state=35):
    pca = PCA(n_components, random_state=random_state) 
    return pca.fit_transform(X)

def norm_data(adata):
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.scale(adata, zero_center=False, max_value=10)
    return adata
    
def hvg (adata):
    sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=3000)
    return adata.var['highly_variable']

def read_adata(path, is_h5ad=False):
    #Todo: add check if Visium, if 'filtered_feature_bc_matrix.h5' or 'data_name_filtered_feature_bc_matrix.h5'
    if is_h5ad:
        adata = sc.read_h5ad(path)
        adata.var_names_make_unique()
        adata.obsm["spatial"]=adata.obsm["spatial"].astype(float)
    else: 
        try:
            adata = sc.read_visium(
        path,
        count_file='filtered_feature_bc_matrix.h5',
        load_images=True
            )

            adata.var_names_make_unique()

            # Ensure spatial coordinates are float
            adata.obsm["spatial"] = (adata.obsm["spatial"].astype(float))

        except Exception as e:

            print(
            f"Warnining scanpy read_visium failed with:\n{e}\n"
            "Falling back to parquet reader..."
            )
            try:
                adata = read_visium_parquet(path)
                adata.var_names_make_unique()

                # Ensure spatial coordinates are float
                adata.obsm["spatial"] = (adata.obsm["spatial"].astype(float))
            except Exception as e:
               import sys
               sys.exit(
                f"Error read visium parquet failed with Error:\n{e}\n"
                )

    return adata



def preprocess(adata, n_components=20,random_seed=35):

    if 'highly_variable' not in adata.var:
        hvg_genes=hvg(adata)
        adata=norm_data(adata)
        adata= adata[:,hvg_genes]
    if 'X_pca' not in adata.obsm:
        if issparse(adata.X):
            data=pca ( adata.X.toarray(), n_components= n_components,random_state=random_seed) 
        else:
            data=pca ( adata.X, n_components= n_components,random_state=random_seed) #if data already a matrix

        adata.obsm["X_pca"]=data
    return adata



import scanpy as sc
import pandas as pd
import numpy as np
import json
from pathlib import Path
from PIL import Image


def read_visium_parquet(
    path,
    count_file="filtered_feature_bc_matrix.h5",
    library_id="visium",
    load_images=True,
):
    """
    Read Visium dataset when spatial coordinates are stored
    in tissue_positions.parquet instead of tissue_positions.csv.

    Parameters
    ----------
    path : str or Path
        Path to Visium sample directory.
    count_file : str
        10x filtered_feature_bc_matrix H5 file.
    library_id : str
        Library identifier stored in adata.uns["spatial"].
    load_images : bool
        Whether to load Visium tissue images.

    Returns
    -------
    adata : AnnData
    """

    path = Path(path)

    # -------------------------
    # Read counts
    # -------------------------
    adata = sc.read_10x_h5(path / count_file)
    adata.var_names_make_unique()

    # -------------------------
    # Read parquet coordinates
    # -------------------------
    tp = pd.read_parquet(
        path / "spatial" / "tissue_positions.parquet"
    )

    tp = tp.set_index("barcode")

    # Match barcode order
    tp = tp.loc[adata.obs_names]

    # Spatial coordinates
    adata.obsm["spatial"] = tp[
        ["pxl_col_in_fullres", "pxl_row_in_fullres"]
    ].to_numpy()

    # Optional metadata
    for col in ["in_tissue", "array_row", "array_col"]:
        if col in tp.columns:
            adata.obs[col] = tp[col].values

    # -------------------------
    # Load images + scalefactors
    # -------------------------
    if load_images:

        spatial_path = path / "spatial"

        # Scalefactors
        with open(
            spatial_path / "scalefactors_json.json"
        ) as f:
            scalefactors = json.load(f)

        images = {}

        hires_path = spatial_path / "tissue_hires_image.png"
        lowres_path = spatial_path / "tissue_lowres_image.png"

        if hires_path.exists():
            images["hires"] = np.array(Image.open(hires_path))

        if lowres_path.exists():
            images["lowres"] = np.array(Image.open(lowres_path))

        adata.uns["spatial"] = {
            library_id: {
                "images": images,
                "scalefactors": scalefactors,
                "metadata": {}
            }
        }

    return adata