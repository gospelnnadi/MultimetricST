

import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.decomposition import PCA
from Evaluate.evaluate import evaluate_cluster
from Evaluate.utils import *
import anndata as ad




import sys, os
ROOT = os.getcwd()
print("The current working dir is, ", ROOT)



def get_adata_from_path(path, is_h5ad):
        import os

        path = args.data_path.rstrip("/")  # remove trailing slash if present

        dir_path = os.path.dirname(path) + "/"  # keep the trailing slash
        base = os.path.basename(path)

        if base.endswith(".h5ad"):
            name = base.rsplit(".", 1)[0]
        else:
            name = base

        print("dir_path:", dir_path)
        print("name:", name)

        adata_raw = read_adata(path, is_h5ad=is_h5ad)
        return adata_raw, dir_path, name
        

def get_adata_from_exp_spatial_path(mode, expression_path, spatial_path):
    if mode==3:   
        spatial_matrix=load_from_file_csv_tsv_npy(args,spatial_path)
        adata = ad.AnnData( )
        adata.obsm["spatial"] = spatial_matrix
        adata.obsm["spatial"]=adata.obsm["spatial"].astype(float)
    else:
        exp_matrix=load_from_file_csv_tsv_npy(args,expression_path)
        spatial_matrix=load_from_file_csv_tsv_npy(args,spatial_path)

        adata = ad.AnnData(X=exp_matrix)
        adata.obsm["spatial"] = spatial_matrix
        adata.obsm["spatial"]=adata.obsm["spatial"].astype(float)
        adata.var_names_make_unique()
    return adata




def plot_label(adata, plot_size, key, savepath):
    import matplotlib.pyplot as plt
    import squidpy as sq
    adata1 = adata[adata.obs[key] != '-1'].copy()
    if plot_size:
      sc.pl.spatial(
        adata1,
        color=key,
        spot_size=plot_size,
        show=False   
        )
      
      plt.savefig(f'{savepath}.png', bbox_inches="tight", dpi=300)
      plt.close()
      adata1.uns.pop(f'{key}_colors')
    else: 
      sc.pl.spatial(adata1, color=key,
                     show=False ) 
      #sq.pl.spatial_scatter(adata, color=key,cmap='Paired', save=savepath) 
      plt.savefig(f'{savepath}.png', bbox_inches="tight", dpi=300)
      plt.close()
      adata1.uns.pop(f'{key}_colors')




def validate_arguments(args):
    """
    Validate input arguments based on the selected mode.
    """
    
    # Mode 1: Requires data for clustering and evaluation
    if args.mode == 1:
       if not args.data_path and (not args.expression_path and not args.spatial_path):
            sys.exit("Error: --data_path to adata object or --expression_path and --spatial_path are required for Mode 1")
                
       if args.data_path and (not args.expression_path and not args.spatial_path):
            if not os.path.exists(args.data_path):
                sys.exit(f"Error: Data file {args.data_path} not found")
            if args.is_h5ad and not ".h5ad" in args.data_path:
                sys.exit(f"Error: Data file {args.data_path} is not in h5ad format specified by --is_h5ad 0 if adat is visium folder containing .h5 file")
            
       if not args.data_path and (args.expression_path):
            if not os.path.exists(args.expression_path):
                sys.exit(f"Error: Data file {args.expression_path} not found")
       if not args.data_path and (args.spatial_path):
            if not os.path.exists(args.spatial_path):
                sys.exit(f"Error: Data file {args.spatial_path} not found")
             
       print("Running full pipeline: clustering → evaluation → dashboard visualization")
       return True
    
    # Mode 2: Requires data and cluster labels for evaluation
    elif args.mode == 2:
       if not args.data_path and (not args.expression_path and not args.spatial_path):
        if not args.data_path:
            sys.exit("Error: --data_path or --expression_path and --spatial_path are required for Mode 2")

       if args.data_path and (not args.expression_path and not args.spatial_path):
            if not os.path.exists(args.data_path):
                sys.exit(f"Error: Data file {args.data_path} not found")
            if args.is_h5ad and not ".h5ad" in args.data_path:
                sys.exit(f"Error: Data file {args.data_path} is not in h5ad format specified by --is_h5ad 0 if adat is visium folder containing .h5 file")
            
       if not args.data_path and (args.expression_path):
            if not os.path.exists(args.expression_path):
                sys.exit(f"Error: Data file {args.expression_path} not found")
       if not args.data_path and (args.spatial_path):
            if not os.path.exists(args.spatial_path):
                sys.exit(f"Error: Data file {args.spatial_path} not found")
             
            
        # For Mode 2, we need either cluster_label_path or method_cluster_labels
       if not args.cluster_label_path and not args.method_cluster_label:
            sys.exit("Error: Either --cluster_label_path or --method_cluster_label is required for Mode 2")
            
        
       if args.cluster_label_path and not os.path.exists(args.cluster_label_path):
            sys.exit(f"Error: Cluster label file {args.cluster_label_path} not found")
            
            
       print("Running Mode 2: evaluation → dashboard visualization")
       return True
    
    # Mode 3: Requires data and evaluation scores CSV for visualization
    elif args.mode == 3:
        args.visual_tissue=1
        if not args.data_path or not args.spatial_path:
            args.visual_tissue = 0
            print("Warning: --data_path or --spatial_path not provided. Spatial visualization will be skipped in Mode 3.")
            
        
        # For Mode 3, we need either spatial data for visualization
        elif args.data_path or args.spatial_path:
            if args.data_path:
                if not os.path.exists(args.data_path):
                    sys.exit(f"Error: Data file {args.data_path} not found")
                if args.is_h5ad and not ".h5ad" in args.data_path:
                    sys.exit(f"Error: Data file {args.data_path} is not in h5ad format specified by --is_h5ad 0 if adat is visium folder containing .h5 file")
            
            else: 
                if not os.path.exists(args.spatial_path):
                    sys.exit(f"Error: Spatial matrix {args.spatial_path} not found")

            if not args.cluster_label_path and not args.method_cluster_label:
                sys.exit("Error: Either --cluster_label_path or --method_cluster_label is required for Mode 3")
            
        
            if args.cluster_label_path and not os.path.exists(args.cluster_label_path):
                sys.exit(f"Error: Cluster label file {args.cluster_label_path} not found")
            elif not os.path.exists(f"{ROOT}/multimetricST_outputs/{args.result_filename}"):
                sys.exit(f"Error: Result scores file {ROOT}/multimetricST_outputs/{args.result_filename} not found")
                
        return True
    
    else:
        print(f"Error: Invalid mode {args.mode}. Mode must be 1, 2, or 3.")


def log_print(mode, data_name):
    import os
    import sys
    import logging

    log_dir = f"{ROOT}/logs"
    os.makedirs(log_dir, exist_ok=True)
    log_file = f"{log_dir}/multimetricst_mode{mode}_{data_name}.log"
    if os.path.exists(log_file):
        os.remove(log_file)

    sys.stdout = PrintLogger(log_file)
    sys.stderr = PrintLogger(log_file)

import sys
from datetime import datetime

class PrintLogger:
    def __init__(self, filename):
        self.file = open(filename, "a")
        self.stdout = sys.__stdout__
        self.buffer = ""

    def write(self, data):
        self.stdout.write(data)  # keep normal console behavior

        self.buffer += data

        while "\n" in self.buffer:
            line, self.buffer = self.buffer.split("\n", 1)

            # Skip empty or progress-bar control lines
            if line.strip() == "" or "\r" in line:
                self.file.write(line + "\n")
            else:
                timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                self.file.write(f"{timestamp} | {line}\n")

        self.file.flush()

    def flush(self):
        self.stdout.flush()
        self.file.flush()
        
    def isatty(self):
        return self.stdout.isatty()

    def fileno(self):
        return self.stdout.fileno()



import argparse
import os
import sys

def main(args):
    """
    Main function that orchestrates the three modes of operation based on user input.
    """
    log_print(args.mode, args.data_name)
    # Validate input parameters based on mode
    if not validate_arguments(args):
        return
    if args.data_path:
        adata_raw, dir_path, data_name = get_adata_from_path(args.data_path, args.is_h5ad)
    else:
        adata_raw=get_adata_from_exp_spatial_path(args.mode , args.expression_path, args.spatial_path)
    data_name=args.data_name
    print("adata_raw dim", adata_raw.shape)
    if os.path.exists(args.ground_truth):
        ground_truth=load_from_file_csv_tsv_npy(args,args.ground_truth, args.ground_truth_col_name)
        print("ground_truth dim", ground_truth.shape)
        # Handle different possible data types (NumPy array or pandas DataFrame/Series)
        if isinstance(ground_truth, (pd.DataFrame, pd.Series)):

            ground_truth =ground_truth.replace(["NA", "nan", ""], "-1")
            ground_truth = np.asarray(ground_truth.replace(np.nan, "-1")).ravel()
            gt = pd.Series(ground_truth)
            mask = gt != '-1'   # keep -1 untouched
            gt_recoded = gt.copy()
            gt_recoded[mask] = (
            gt[mask]
           .astype("category")
           .cat.codes + 1
           ).astype(str)
            ground_truth = gt_recoded.to_numpy()
        elif isinstance(ground_truth, np.ndarray):
            ground_truth = pd.Series(ground_truth).replace(["NA", "nan", ""], np.nan).to_numpy()
            ground_truth = np.where(np.isnan(ground_truth), "-1", ground_truth).ravel()
            gt = pd.Series(ground_truth)
            mask = gt != '-1'   # keep -1 untouched
            gt_recoded = gt.copy()
            gt_recoded[mask] = (
            gt[mask]
           .astype("category")
           .cat.codes + 1
           ).astype(str)
            ground_truth = gt_recoded.to_numpy()
        else:
            # Fallback: if it's a list or other type
            ground_truth = np.array(ground_truth)
            ground_truth = pd.Series(ground_truth).replace(["NA", "nan", ""], np.nan).to_numpy()
            ground_truth = np.where(np.isnan(ground_truth), "-1", ground_truth).ravel()
            gt = pd.Series(ground_truth)
            mask = gt != '-1'   # keep -1 untouched
            gt_recoded = gt.copy()
            gt_recoded[mask] = (
            gt[mask]
           .astype("category")
           .cat.codes + 1
           ).astype(str)
            ground_truth = gt_recoded.to_numpy()
        
        print("ground_truth unique", np.unique(ground_truth[ground_truth != '-1']))
        adata_raw.obs["ground_truth"] = ground_truth.astype(str)
    elif args.ground_truth in adata_raw.obs: 
        ground_truth = adata_raw.obs[args.ground_truth]
        if isinstance(ground_truth, (pd.DataFrame, pd.Series)):
            ground_truth =ground_truth.replace(["NA", "nan", ""], "-1")
            ground_truth = np.asarray(ground_truth.replace(np.nan, "-1")).ravel()
            gt = pd.Series(ground_truth)
            mask = gt != '-1'   # keep -1 untouched
            gt_recoded = gt.copy()
            gt_recoded[mask] = (
            gt[mask]
           .astype("category")
           .cat.codes + 1
           ).astype(str)
            ground_truth = gt_recoded.to_numpy()
        elif isinstance(ground_truth, np.ndarray):
            ground_truth = pd.Series(ground_truth).replace(["NA", "nan", ""], np.nan).to_numpy()
            ground_truth = np.where(np.isnan(ground_truth), "-1", ground_truth).ravel()
            gt = pd.Series(ground_truth)
            mask = gt != '-1'  # keep -1 untouched
            gt_recoded = gt.copy()
            gt_recoded[mask] = (
            gt[mask]
           .astype("category")
           .cat.codes + 1
           ).astype(str)
            ground_truth = gt_recoded.to_numpy()
        else:
            # Fallback: if it's a list or other type
            ground_truth = np.array(ground_truth)
            ground_truth = pd.Series(ground_truth).replace(["NA", "nan", ""], np.nan).to_numpy()
            ground_truth = np.where(np.isnan(ground_truth), "-1", ground_truth).ravel()
            gt = pd.Series(ground_truth)
            mask = gt != '-1'  # keep -1 untouched
            gt_recoded = gt.copy()
            gt_recoded[mask] = (
            gt[mask]
           .astype("category")
           .cat.codes + 1
           ).astype(str)
            ground_truth = gt_recoded.to_numpy()
        adata_raw.obs['ground_truth']=ground_truth

        adata_raw.obs["ground_truth"] = ground_truth.astype(str)
        print("ground_truth dim", ground_truth.shape)
        print("ground_truth unique", np.unique(ground_truth[ground_truth != '-1']))
    else:
        print("Warning: Ground truth annotation not found in provided path or adata.obs. External metrics will be set to 0.")
        args.external_metrics=False
        if "ground_truth" in adata_raw.obs.columns:
            del adata_raw.obs["ground_truth"]
      
        

    # Mode 1: Full pipeline - run clustering, evaluation, and dashboard
    if args.mode == 1:
        print("Running Mode 1: Full pipeline")
        run_full_pipeline(args,adata_raw,data_name,subset_methods=args.subset_methods,data_type=args.data_type,n_clusters=args.n_clusters)
        print("Full pipeline completed.")
    
    # Mode 2: Evaluation and visualization only
    elif args.mode == 2:
        print("Running Mode 2: Evaluation and visualization")
        run_evaluation_and_visualization(args,adata_raw,data_name, data_type=args.data_type)
        print("Evaluation and visualization completed.")
    
    # Mode 3: Visualization only from precomputed scores
    elif args.mode == 3:
        print("Running Mode 3: Visualization from precomputed scores")
        run_visualization_only(args,adata_raw,data_name)
        print("Visualization completed.")



def run_full_pipeline(args,adata_raw,data_name,subset_methods=None,data_type='Visium',n_clusters=7,n_components=20, random_seed=35,):
    """
    Execute Mode 1: Run clustering, evaluation, and visualization.
    """
    try:
        # Import and run clustering module
       
        from Spatial_Clustering_Methods.clustering import run_clustering_pipeline
        print("Running clustering pipeline...")
        
        # Run clustering and get cluster labels
        adata, comp_cost = run_clustering_pipeline(adata_raw.copy(),data_name,subset_methods=subset_methods,data_type=args.data_type,n_clusters=n_clusters)
        
        # Import and run evaluation module
        from Evaluate.evaluate import evaluate_cluster
        print("Running evaluation pipeline...")
        # Run evaluation and get scores
        adata= preprocess(adata, n_components=n_components, random_seed=random_seed)

        plot_savepath=f"{ROOT}/multimetricST_outputs/figures/{data_name}/"
        os.makedirs(plot_savepath, exist_ok=True)
        if 'ground_truth' in adata.obs:
            ground_truth =adata.obs['ground_truth'] 
            print("ground truth was detcted in anndata and will be used for evaluation")
            method='ground_truth'
            print(f"Saving {method} cluster plot to {plot_savepath}")
            plot_label(adata, plot_size=args.plot_size, key=method, savepath=f"{plot_savepath}/{method}" )
            print(f"Saved {method} cluster plot")
        else: 
            ground_truth=None
            print("ground truth was detcted in anndata and will be used for evaluation")
        
        # PCA matrix
        pca_matrix = adata.obsm["X_pca"]
        if data_type=="Visium":
            is_visium=True
        else:
            is_visium=False
        results=[]

        if args.external_metrics:
            print("Evaluate external metrics is True...")
            if ground_truth is None or len(ground_truth) == 0:
                print("Warning: Ground truth annotation is missing or empty. External metrics will be set to 0.")   
            else:
                print(f"Ground truth annotation provided with {len(np.unique(ground_truth))} unique clusters.")
                
        if args.internal_metrics:
            print("Evaluate internal metrics is True...")
        elif args.internal_metrics == False:
            print("Evaluate internal metrics is False. Only external metrics will be computed if ground truth is available.")
        
        for m in comp_cost:
            method=m['method']
            # Clusters from adata
            print(f"\n\nEvaluating clustering results...{method}\n\n")

            if adata[adata.obs[method] == '-1'].shape[0] > 0:
                adata1 = adata[adata.obs[method] != '-1'].copy()
                adata1 = preprocess(adata1, n_components=n_components, random_seed=random_seed)
                pca_matrix1 = adata1.obsm["X_pca"]
 
            else:   
                adata1= adata.copy()
                pca_matrix1 = pca_matrix
                
            pred = adata1.obs[method]
            if 'ground_truth' in adata1.obs:
                    ground_truth =adata1.obs['ground_truth'] 
        
            scores=evaluate_cluster( adata1,pred, ground_truth, pca_matrix1, is_visium=is_visium,verbose=True,decimal=4, external_metrics=args.external_metrics, internal_metrics=args.internal_metrics)
            
            # Merge metrics and computational cost into one record
            m.update(scores)

            # Append to list
            results.append(m)
            print(f"\n\nSaving {method} cluster plot to {plot_savepath}\n\n")
            plot_label(adata, plot_size=args.plot_size, key=method, savepath=f"{plot_savepath}/{method}" )
            print(f"Saved {method} cluster plot")
        # Convert to DataFrame
        results_df = pd.DataFrame(results) #.round(4)
        result_savepath=f"{ROOT}/multimetricST_outputs/{args.result_filename}"
        results_df.to_csv(result_savepath,  index=False)
        for col in adata.obs.columns:
            if adata.obs[col].dtype == 'object':
                adata.obs[col] = adata.obs[col].astype(str)
        os.makedirs("Data/Preprocessed", exist_ok=True)
        adata.write_h5ad(f"Data/Preprocessed/{data_name}.h5ad")
        print("Clustering and evaluation results saved.")

        print("Running dashboard visualization...")
        # Run dashboard visualization
        run_dashboard(result_savepath, plot_savepath, port=args.dashboard_port )
 

    except Exception as e:
        print(f"Error in full pipeline execution: {str(e)}")
        raise
        

def run_evaluation_and_visualization(args,adata_raw,data_name,data_type='Visium',n_components=20, random_seed=35):
    """
    #Execute Mode 2: Run evaluation and visualization only.
    """ 
    try:
        # Load cluster labels
        print(" show methods_cluster_label ", args.method_cluster_label)
        print(" show methods_cluster_label type ", type(args.method_cluster_label))

        cluster_labels = None
        if args.cluster_label_path and os.path.exists(args.cluster_label_path):
            cluster_labels = load_from_file_csv_tsv_npy(args,args.cluster_label_path, args.cluster_label_col_names)
            args.method_cluster_label = cluster_labels.columns.tolist()
            for col in cluster_labels.columns:
                adata_raw.obs[col] = cluster_labels[col]
        elif all(label in adata_raw.obs for label in args.method_cluster_label):
            print("Using method_cluster_label from anndata.obs for evaluation")
        else:
            sys.exit("Error: No cluster labels provided for evaluation in Mode 2")
        
        from Evaluate.evaluate import evaluate_cluster
        print("Running evaluation pipeline...")
        # Run evaluation and get scores
        adata_raw = preprocess(adata_raw, n_components=n_components, random_seed=random_seed)

        plot_savepath=f"{ROOT}/multimetricST_outputs/figures/{data_name}/"
        os.makedirs(plot_savepath, exist_ok=True)
        if 'ground_truth' in adata_raw.obs:
            ground_truth =adata_raw.obs['ground_truth'] 
            print("ground truth was detcted in anndata and will be used for evaluation")
            method='ground_truth'
            print(f"Saving {method} cluster plot to {plot_savepath}")
            plot_label(adata_raw, plot_size=args.plot_size, key=method, savepath=f"{plot_savepath}/{method}" )
            print(f"Saved {method} cluster plot")
        else: 
            ground_truth=None
            print("ground truth was detcted in anndata and will be used for evaluation") 
        
        # PCA matrix
        pca_matrix = adata_raw.obsm["X_pca"]
        if data_type=="Visium":
            is_visium=True
        else:
            is_visium=False
        results=[]
        if args.external_metrics:
            print("Evaluate external metrics is True...")
            if ground_truth is None or len(ground_truth) == 0:
                print("Warning: Ground truth annotation is missing or empty. External metrics will be set to 0.")   
            else:
                print(f"Ground truth annotation provided with {len(np.unique(ground_truth))} unique clusters.")
                
        if args.internal_metrics:
            print("Evaluate internal metrics is True...")
        elif args.internal_metrics == False:
            print("Evaluate internal metrics is False. Only external metrics will be computed if ground truth is available.")
  
        for method in args.method_cluster_label:
            if method not in adata_raw.obs:
                print(f"Warning: Method {method} not found in adata.obs. Skipping evaluation and plot.")
                continue
            else: 
                print(f"\n\nEvaluating clustering results...{method}\n\n")
            # Clusters from adata
            if adata_raw[adata_raw.obs[method] == '-1'].shape[0] > 0:
                adata1 = adata_raw[adata_raw.obs[method] != '-1'].copy()
                adata1 = preprocess(adata1, n_components=n_components, random_seed=random_seed)
                pca_matrix1 = adata1.obsm["X_pca"]
     
            else:   
                adata1= adata_raw.copy()
                pca_matrix1 = pca_matrix
            pred = adata1.obs[method]
            if 'ground_truth' in adata1.obs:
                    ground_truth =adata1.obs['ground_truth'] 

            scores=evaluate_cluster( adata1,pred, ground_truth, pca_matrix1, is_visium=is_visium,verbose=True,decimal=4, external_metrics=args.external_metrics, internal_metrics=args.internal_metrics)
            
                # Merge metrics and computational cost into one record
            m={'method':method}
            m.update(scores)

            # Append to list
            results.append(m)
            print(f"\n\nSaving {method} cluster plot to {plot_savepath}\n\n")
            plot_label(adata_raw, plot_size=args.plot_size, key=method, savepath=f"{plot_savepath}/{method}" )
            print(f"Saved {method} cluster plot")

        # Convert to DataFrame
        results_df = pd.DataFrame(results).round(4)
        result_savepath=f"{ROOT}/multimetricST_outputs/{args.result_filename}"
        results_df.to_csv(result_savepath,  index=False)

        print("Evaluation results saved.")
        print("Running dashboard visualization...")
        # Run dashboard visualization
        run_dashboard(result_savepath, plot_savepath, port=args.dashboard_port )
        
        
    except Exception as e:
        print(f"Error in evaluation and visualization execution: {str(e)}")
        raise 

def run_visualization_only(args,adata_raw,data_name):
    """
    #Execute Mode 3: Visualization only from precomputed scores.
    """
    try:
        # For Mode 3, we need to load precomputed evaluation scores
        # This would typically be from a CSV file specified by the user
        # For now, we'll assume the scores are generated in previous steps
        # or provided through additional arguments
        
        # Load cluster labels if provided for spatial visualization
        file_path=f"{ROOT}/multimetricST_outputs/{args.result_filename}"
        cluster_labels = None
        if args.visual_tissue and args.cluster_label_path and os.path.exists(args.cluster_label_path):
            cluster_labels = load_from_file_csv_tsv_npy(args,args.cluster_label_path, args.cluster_label_col_names)
            args.method_cluster_label = list(set(args.method_cluster_label) | set(cluster_labels.columns.tolist()))
            for col in cluster_labels.columns:
                adata_raw.obs[col] = cluster_labels[col]
        
        elif file_path and os.path.exists(file_path):
            args.method_cluster_label  = load_from_file_csv_tsv_npy(args,file_path, None)
            print("Using method_cluster_label from precomputed results for visualization")
            args.method_cluster_label = args.method_cluster_label['method'].tolist()
            print(" show methods_cluster_label ", args.method_cluster_label)
        

        plot_savepath=f"{ROOT}/multimetricST_outputs/figures/{data_name}/"
        os.makedirs(plot_savepath, exist_ok=True)
        if 'ground_truth' in adata_raw.obs:
        
            method='ground_truth'
            print(f"Saving {method} cluster plot to {plot_savepath}")
            plot_label(adata_raw, plot_size=args.plot_size, key=method, savepath=f"{plot_savepath}/{method}" )
            print(f"Saved {method} cluster plot")
        
        if args.method_cluster_label and adata_raw is not None:
            for method in args.method_cluster_label:
                print(f"Saving {method} cluster plot to {plot_savepath}")
                if method in adata_raw.obs:
                    plot_label(adata_raw, plot_size=args.plot_size, key=method, savepath=f"{plot_savepath}/{method}"  )
                    print(f"Saved {method} cluster plot")
                else:
                    print(f"Warning: Method {method} not found in adata.obs. Skipping plot.")
      

        # Convert to DataFrame
        
        result_savepath=f"{ROOT}/multimetricST_outputs/{args.result_filename}"
        print("Running dashboard visualization...")
        # Run dashboard visualization
        run_dashboard(result_savepath,plot_savepath, port=args.dashboard_port )
        

    except Exception as e:
        print(f"Error in visualization-only execution: {str(e)}")
        raise 

def load_from_file_csv_tsv_npy(args,path, col_name="annotation"):
    """
    Load from file (CSV or TSV or npy format).
    """
    import numpy as np
    import pandas as pd
    is_vector=True
    if (path == args.expression_path or path==args.spatial_path) or col_name==None:
        is_vector=False
    
    if  path.endswith('.npy'):
        return np.load( path)
    elif  path.endswith('.csv'):
        df = pd.read_csv( path)
        # Assume first column contains labels
        if is_vector:
            return df[[c for c in col_name if c in df.columns]].astype(str) #df[col_name].astype(str)
        else:
            return df
    elif path.endswith('.tsv'):
        df= pd.read_csv(path, sep='\t')
        if is_vector:
            return df[[c for c in col_name if c in df.columns]].astype(str)
        else:
            return df
    else:
        raise ValueError(f"Unsupported file format for cluster labels: { path}")

def release_port(PORT=5009):
    import psutil

    for conn in psutil.net_connections(kind='inet'):
      if conn.laddr.port == PORT and conn.pid is not None:
        try:
            proc = psutil.Process(conn.pid)
            print(f"Killing PID {conn.pid} ({proc.name()}) using port {PORT}")
            proc.kill()
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            pass
    print(f"Port {PORT} is now free.")
import socket

def is_port_free(port):
    s = socket.socket()
    try:
        s.bind(("0.0.0.0", port))
        return True
    except:
        return False
    finally:
        s.close()



def run_dashboard(result_savepath,plot_savepath=None, port=5009):
    """
    Run the dashboard visualization module.
    """
    try:
        from Visualize_Scores.dashboard import create_dashboard
       
        release_port(port)  # Ensure port is free before launching
        print("Port free?", is_port_free(port))
        print("Launching dashboard...")
        # Create and launch dashboard
        dashboard = create_dashboard(result_savepath, plot_savepath)
        
        dashboard.show(address="0.0.0.0",port=port, websocket_origin="*",
        session_token_expiration=3600*24  # 24 hour expiration for session tokens 
        )

        print(f"Launching server at http://localhost:{port}")
        
    except Exception as e:
        print(f"Error launching dashboard: {str(e)}")
        raise

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="MultimetricST: A comprehensive framework for spatial transcriptomics clustering evaluation",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    
    # Core execution parameters
    parser.add_argument("--mode", type=int, default=1, choices=[1, 2, 3],
                       help="Execution mode: (1) full pipeline, (2) evaluation and visualization only, (3) visualization only from precomputed scores")
    parser.add_argument("--subset_methods",nargs='+', default=["all"], 
                       help="provide a list of methods to subset the model execution, evaluation and visualization. E.g., GIST GraphST SEDR. Default is all methods")
    
    
    # Data input parameters
    parser.add_argument("--data_path", type=str, default="",
                       help="Full path to the data object. Can be: Data/DLPFC/151673 (folder) or Data/MK.h5ad (file)")
    parser.add_argument("--data_name", type=str, default="DLPFC_151673",
                       help="Dataset Name. Can be: 151673 (DLPFC slice) or Mouse_Kidney etc. ")
    
    parser.add_argument("--data_type", type=str, default='Visium',
                       help="Data technology Visium, Stereo, Slide, Xenium. Default is Visium")
    parser.add_argument("--is_h5ad", type=int, default=0, choices=[0, 1],
                       help="Data is in h5ad format (1) or not (0). Default is 0 for Visium data folder")
    
    # Alternative data input (for non-h5ad data)
    parser.add_argument("--expression_path", type=str, default="",
                       help="Full path to the .npy file containing expression matrix (required for non-h5ad data)")
    parser.add_argument("--spatial_path", type=str, default="",
                       help="Full path to the .npy file containing spatial coordinates (required for non-h5ad data)")
    # Groundtruth annotation 
    parser.add_argument("--ground_truth", type=str, default="ground_truth",
                       help="if the groundtruth annotation is in the anndata specify the anndata.obs key ground_truth. Otherwise, provide the path to the annotation tsv or csv file. " \
                       "If file or obs key not exist Annotation-Dependent evaluation metrics will not be used")
    parser.add_argument("--ground_truth_col_name",nargs='+', default=["layer_guess"], 
                       help="if the path to the annotation tsv or csv file is provided, indicate the column indice where the ground-truth annotation is found. For DLPFC data, it is layer_guess")
    
    # Cluster labels input (for modes 2 and 3)
    parser.add_argument("--cluster_label_path", type=str, default="",
                       help="Full path to .CSV or .npy file containing cluster labels for each spot/cell")
    parser.add_argument("--cluster_label_col_names", nargs='+', default=[],
                       help="if the path to the cluster label tsv or csv file is provided, indicate the column indice where the cluster label is found. E.g., GIST GraphST SEDR")
    
    parser.add_argument("--n_clusters", type=int, default=7,
                       help="Number of clusters")
    # Multiple method cluster labels
    parser.add_argument("--method_cluster_label", nargs='+', default=[],
                       help="Specify the anndata.obs keys of each method cluster label to be evaluated (e.g., GIST GraphST SEDR)")
    
    parser.add_argument("--seed", type=int, default=35,
                       help="Random seed for reproducibility")
    
    # Visualization parameters
    parser.add_argument("--plot_size", type=int, default=0,
                       help="Tissue spatial plot size. Set to 0 when histology image is available in adata")
    parser.add_argument("--result_filename", type=str, default="clustering_results.csv", 
                       help="CSV file name of precomputed clustering results (required for Mode 3), optional for Mode 1 and 2. . The defualt is clustering_results.csv. The file should be saved placed in MultimetricST/multimetricST_outputs/ for Mode 3. The computed results with Mode 1 and 2 will be saved in the same folder.")
    parser.add_argument("--dashboard_port", type=int, default=8008, 
                       help="Specify the port of which to launch the dashboard. default is 8008. Make sure the port is free or change to a different port if needed.")
  
    parser.add_argument("--external_metrics", type=int, default=1, choices=[0, 1],
                       help="Whether to compute external evaluation metrics that require ground truth annotation (1) or not")
    parser.add_argument("--internal_metrics", type=int, default=1, choices=[0, 1],
                       help="Whether to compute internal (annotation independent) evaluation metrics (1) or not")

    args = parser.parse_args()
    main(args)


