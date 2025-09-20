
import numpy as np
from sklearn import metrics
from Silhouette_Spatial_Score.silhouette_spatial import silhouette_spatial_score
from .utils import *  # Assumes CHAOS, PAS, ASW functions are defined here
from scipy.sparse import issparse 
import matplotlib.pyplot as plt
import seaborn as sns



















def evaluate_cluster(adata, pred, ground_truth=None,  pca_matrix=None, is_visium=False,verbose=True,decimal=4):
  
    """
    Evaluate clustering performance using both traditional and spatially-aware metrics.

    This function computes a range of clustering evaluation scores given predicted labels and (optionally) ground truth.
    It supports both standard metrics (e.g., ARI, AMI, silhouette) and spatial transcriptomics-aware metrics such as 
    spatial silhouette score (SSS), CHAOS, PAS, and ASW. Useful for benchmarking clustering quality in spatial omics.

    Parameters
    ---------- 
    adata : AnnData
        AnnData object containing spatial and expression data. Must contain `obsm['spatial']` and optionally 
        `uns['average_penalty']` if spatial penalty score is used.

 
    pred : array-like
        Predicted cluster labels to be evaluated.

    ground_truth : array-like (default: None) if groundth truth annotation is un available
        Ground truth cluster labels. Required for ARI, AMI, Purity, Homogeneity, Completeness, and V-measure.
    
  
    pca_matrix : np.ndarray or None, optional (default: None)
        PCA-reduced expression matrix to be used for distance-based metrics. If None, `adata.X` is used.

    is_visium : bool, optional (default: False)
        Whether the dataset is from 10X Visium platform. Affects how spatial silhouette is computed.

    verbose : bool, optional (default: True)
        Whether to print metric values to standard output.

    decimal : int, optional (default: 4)
        Number of decimal places to round metric values.

    Returns
    -------
    metrics_dict : dict
        Dictionary containing the following keys and their corresponding scores:
            - "ARI"              : Adjusted Rand Index
            - "AMI"              : Adjusted Mutual Information
            - "Purity"           : Cluster purity
            - "Homogeneity"      : Homogeneity score
            - "Completeness"     : Completeness score
            - "V-Measure"        : Harmonic mean of homogeneity and completeness
            - "SSS"              : Spatial silhouette score (Silhouette_Spatial_Score)
            - "SSS-Penalty"      : Precomputed spatial penalty (from adata.uns['average_penalty'])
            - "Silhouette"       : Standard silhouette score (cosine distance)
            - "Davies-Bouldin"   : Davies-Bouldin index (lower is better)
            - "CHAOS"            : Cluster spatial compactness (%)
            - "PAS"              : Preservation of spatial autocorrelation
            - "ASW"              : Average Spatial Width

    Notes
    -----
    - CHAOS, PAS, and ASW are custom spatial metrics typically defined in `.utils`.
    - If fewer than 2 predicted clusters are provided, all metrics default to 0.
    - Assumes `adata.obsm['spatial']` is present and contains 2D spatial coordinates for each spot.
    - Assumes `adata.uns['average_penalty']` is computed by SSS.

    Example
    -------
    >>> scores = evaluate_cluster(adata=adata, pred=y_pred, ground_truth=y_true, pca_matrix=pca, is_visium=False)
    >>> print(scores["ARI"], scores["SSS"], scores["CHAOS"])
    """

    # If no PCA matrix is provided, use adata.X (convert from sparse if necessary)
    if pca_matrix is None:
        if issparse(adata.X):
            pca_matrix = adata.X.toarray()
        else:
            pca_matrix = adata.X

    # Initialize scores
    if ground_truth is not None and len(ground_truth):
        # Adjusted Rand Index (ARI)
        ARI = metrics.adjusted_rand_score(ground_truth, pred)
        ARI = np.round(ARI, decimal)

        # Adjusted Mutual Information (AMI)
        AMI = metrics.adjusted_mutual_info_score(ground_truth, pred)
        AMI = np.round(AMI, decimal)

        # Purity Score: sum of max over rows in contingency matrix / total
        purity = metrics.cluster.contingency_matrix(ground_truth, pred).max(axis=1).sum() / len(ground_truth)
        purity = np.round(purity, decimal)

        # Homogeneity Score
        homogeneity = metrics.homogeneity_score(ground_truth, pred)
        homogeneity = np.round(homogeneity, decimal)

        # Completeness Score
        completeness = metrics.completeness_score(ground_truth, pred)
        completeness = np.round(completeness, decimal)

        # V-measure (harmonic mean of homogeneity and completeness)
        v_measure = metrics.v_measure_score(ground_truth, pred)
        v_measure = np.round(v_measure, decimal)
    else:
        # Fallback to 0 if no valid ground truth
        ARI, AMI, purity, homogeneity, completeness, v_measure = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0

    # Compute clustering metrics if at least 2 clusters are present
    if len(np.unique(pred)) > 1:
        # Spatial silhouette score (custom metric considering spatial info)
        silhouette_spatial = silhouette_spatial_score(pca_matrix, pred, adata, metric="cosine", is_visium=is_visium)
        silhouette_spatial = np.round(silhouette_spatial, decimal)

        # Penalty based on spatial smoothness or inconsistency (precomputed)
        penalty = np.round(adata.uns['average_penalty'], decimal)

        # Standard silhouette score with cosine distance
        silhouette = metrics.silhouette_score(pca_matrix, pred, metric='cosine')
        silhouette = np.round(silhouette, decimal)

        # Davies-Bouldin index (lower is better)
        davies_bouldin = metrics.davies_bouldin_score(pca_matrix, pred)
        davies_bouldin = np.round(davies_bouldin, decimal)

        # CHAOS: cluster homogeneity metric in spatial context
        chaos = compute_CHAOS(pred, adata.obsm['spatial'])
        chaos = np.round(chaos * 100 , decimal) # scale to percentage

        # PAS: spatial autocorrelation preservation metric
        pas = compute_PAS(pred, adata.obsm['spatial'])
        pas = np.round(pas, decimal)

        # ASW: average spatial within-cluster distance
        ASW = compute_ASW(pred, adata.obsm['spatial'])
        ASW = np.round(ASW, decimal)

    else:
        # Fallback values for single-cluster predictions
        silhouette_spatial, silhouette, davies_bouldin, chaos, pas, ASW = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
        print("Cluster size is less than 2")

    # Combine all results into a list
    eval_result = [
        ARI, AMI, purity, homogeneity, completeness, v_measure,
        silhouette_spatial, penalty, silhouette, davies_bouldin,
        chaos, pas, ASW
    ]

    # Print metrics if verbose mode is on
    if verbose:
        print('ARI: ', ARI)
        print("AMI: ", AMI)
        print("Purity Score: ", purity)
        print("Homogeneity Score: ", homogeneity)
        print("Completeness Score: ", completeness)
        print("V-Measure Score: ", v_measure)
        print("Silhouette Spatial Score: ", silhouette_spatial)
        print("SSS average_penalty:", penalty)
        print("Silhouette Score: ", silhouette)
        print("Davies Bouldin Index: ", davies_bouldin)
        print("CHAOS: ", chaos)
        print("PAS: ", pas)
        print("ASW: ", ASW)

    # Construct and return a dictionary of metrics
    metrics_dict = {
        "ARI": eval_result[0],
        "AMI": eval_result[1],
        "Purity": eval_result[2],
        "Homogeneity": eval_result[3],
        "Completeness": eval_result[4],
        "V-Measure": eval_result[5],
        "Silhouette-Spatial": eval_result[6],
        "SSS-Penalty": eval_result[7],
        "Silhouette": eval_result[8],
        "Davies-Bouldin": eval_result[9],
        "CHAOS": eval_result[10],
        "PAS": eval_result[11],
        "ASW": eval_result[12]
    }

    return metrics_dict




