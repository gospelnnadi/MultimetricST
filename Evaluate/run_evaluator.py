

from sklearn.decomposition import PCA
from Evaluate.evaluate import evaluate_cluster
from Evaluate.utils import *

import warnings
import anndata as ad


import warnings
import anndata as ad

def evaluate_adata(
    model_adata,
    ground_truth=None,
    pred=None,
    path_raw_adata="",
    dataname_raw_adata="",
    is_h5ad=False,
    n_components=20,
    random_seed=35,
    is_visium=True,
):
    """
    Evaluate clustering performance given an existing AnnData object 
    (with optional raw dataset loading).

    Parameters
    ----------
    model_adata : AnnData
        AnnData object containing expression/spatial data and possibly predictions.
    ground_truth : array-like, optional
        Ground-truth labels. If not provided, a warning is issued.
    pred : array-like, optional
        Predicted cluster labels. If None, will attempt to use model_adata.obs["cluster"].
    path_raw_adata : str, optional
        Path to the folder containing the raw AnnData file.
    dataname_raw_adata : str, optional
        File name (without extension) of the raw AnnData file.
    is_h5ad : bool, default=False
        Whether the raw AnnData file is stored in `.h5ad` format.
    n_components : int, default=20
        Number of PCA components to compute in preprocessing.
    random_seed : int, default=35
        Random seed for reproducibility.
    is_visium : bool, default=True
        Whether the data is Visium-type spatial transcriptomics.

    Returns
    -------
    scores : dict
        Dictionary of evaluation metrics computed by `evaluate_cluster`.
    """

    # --- Decide whether to read raw AnnData or use model_adata ---
    path_data = True if path_raw_adata != "" and dataname_raw_adata != "" else False
    if not path_data:
        print("path empty using model_adata directly")
        if not isinstance(model_adata, ad.AnnData):
            raise ValueError("model_adata must be an AnnData object")
        adata_raw = model_adata
    else:
        adata_raw = read_adata(
            f"{path_raw_adata}/{dataname_raw_adata}", is_h5ad=is_h5ad
        )

    # --- Preprocess dataset (e.g., PCA embedding) ---
    adata_raw = preprocess(adata_raw, n_components=n_components, random_seed=random_seed)

    # --- Handle predictions ---
    if pred is None:
        if "cluster" in model_adata.obs:
            adata_raw.obs["cluster"] = model_adata.obs["cluster"]
            pred = "cluster"
        else:
            raise ValueError(
                "No predictions provided and 'cluster' not found in model_adata.obs"
            )
    else:
        adata_raw.obs["cluster"] = pred
        pred = "cluster"

    # --- Handle ground truth ---
    if ground_truth is None:
        warnings.warn("No ground_truth provided", UserWarning)
    else:
        adata_raw.obs["ground_truth"] = ground_truth
        ground_truth = "ground_truth"

    # --- Align indices if raw data was loaded separately ---
    if path_data and model_adata is not None:
        common_idx = adata_raw.obs_names.intersection(model_adata.obs_names)
        adata_raw = adata_raw[common_idx].copy()

    # --- Extract PCA matrix ---
    pca_matrix = adata_raw.obsm["X_pca"]

    # --- Evaluate clustering ---
    scores = evaluate_cluster(
        adata_raw,
        pred,
        ground_truth,
        pca_matrix,
        is_visium=is_visium,
        verbose=True,
        decimal=4,
    )
    return scores


def evaluate_X(
    exp_matrix,
    spatial_matrix,
    pred=None,
    ground_truth=None,
    n_components=20,
    random_seed=35,
    is_visium=True,
):
    """
    Evaluate clustering performance given raw expression and spatial matrices.

    Parameters
    ----------
    exp_matrix : np.ndarray or scipy.sparse
        Expression matrix (cells/spots × genes).
    spatial_matrix : np.ndarray
        Spatial coordinates matrix (cells/spots × 2).
    pred : array-like, optional
        Predicted cluster labels. Required.
    ground_truth : array-like, optional
        Ground-truth labels. If not provided, a warning is issued.
    n_components : int, default=20
        Number of PCA components to compute in preprocessing.
    random_seed : int, default=35
        Random seed for reproducibility.
    is_visium : bool, default=True
        Whether the data is Visium-type spatial transcriptomics.

    Returns
    -------
    scores : dict
        Dictionary of evaluation metrics computed by `evaluate_cluster`.
    """

    # --- Dimension checks ---
    n_exp = exp_matrix.shape[0]
    n_spatial = spatial_matrix.shape[0]
    n_pred = len(pred) if pred is not None else None
    n_gt = len(ground_truth) if ground_truth is not None else None

    if n_exp != n_spatial:
        raise ValueError(
            f"exp_matrix ({n_exp}) and spatial_matrix ({n_spatial}) must have same number of cells/spots"
        )
    if n_pred is not None and n_pred != n_exp:
        raise ValueError(f"pred length ({n_pred}) does not match n_exp ({n_exp})")
    if n_gt is not None and n_gt != n_exp:
        raise ValueError(
            f"ground_truth length ({n_gt}) does not match n_exp ({n_exp})"
        )

    # --- Compose AnnData object from expression + spatial ---
    adata = ad.AnnData(X=exp_matrix)
    adata.obsm["spatial"] = spatial_matrix

    # --- Preprocess (e.g., PCA embedding) ---
    adata_proc = preprocess(adata, n_components=n_components, random_seed=random_seed)

    # --- Predictions ---
    if pred is None:
        raise ValueError("No predictions provided")
    else:
        adata_proc.obs["cluster"] = pred
        pred_key = "cluster"

    # --- Ground truth ---
    if ground_truth is None:
        warnings.warn("No ground_truth provided", UserWarning)
        ground_truth_key = None
    else:
        adata_proc.obs["ground_truth"] = ground_truth
        ground_truth_key = "ground_truth"

    # --- PCA matrix ---
    pca_matrix = adata_proc.obsm["X_pca"]

    # --- Evaluate clustering ---
    scores = evaluate_cluster(
        adata_proc,
        pred_key,
        ground_truth_key,
        pca_matrix,
        is_visium=is_visium,
        verbose=True,
        decimal=4,
    )
    return scores
