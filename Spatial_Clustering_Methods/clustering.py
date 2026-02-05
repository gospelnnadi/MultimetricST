
import subprocess
import os


import sys, os
current_dir = f"{os.getcwd()}/Spatial_Clustering_Methods"
ROOT = f"{os.path.dirname(current_dir)}/Spatial_Clustering_Methods" #os.getcwd() #os.path.dirname(os.path.abspath(__file__))  # project root

# Add all method directories to sys.path
method_dirs = [
    
    "DeepST",
    #"GraphST",
    "GIST",
    "SEDR",
    "SpaGCN",
    "STAGATE_pyG",
    "SCAN-IT",
    "CCST",
    "conST",     
    "SpaceFlow"
]

for d in method_dirs:
    path = os.path.join(current_dir, d)
    if path not in sys.path:
        sys.path.append(path)
sys.path.append(ROOT)







from runDeepST import run as runDeepST
from runGraphST import run as runGraphST
from runGIST import run as runGIST
from runSEDR import run as runSEDR 
from runSpaGCN import run as runSpaGCN
from runSTAGATE import run as runSTAGATE 
from runScanIT import run as runScanIT
from runCCST import run as runCCST
from runConST import run as runConST 
from runSpaceFlow import run as runSpaceFlow 

import pandas as pd
import numpy as np

"""   ("DeepST", runDeepST) ,
    ("CCST", runCCST),
    ("conST", runConST),
    ("GIST", runGIST), 
    ("GraphST", runGraphST), 
    ("SCAN-IT", runScanIT), 
    ("SEDR", runSEDR),
    ("SpaceFlow", runSpaceFlow) ,
    ("SpaGCN", runSpaGCN), 
    ("STAGATE", runSTAGATE)  """
methods = [ 
    ("SpaGCN", runSpaGCN), 
    ("STAGATE", runSTAGATE),
    ("SpaceFlow", runSpaceFlow)
]
comp_cost = []
def run_clustering_pipeline(adata_raw,data_name,data_type='Visium',n_clusters=7,decimal=4):
  for method_name, method_func in methods:
    print(f"\n\nRunning {method_name}...\n\n\n\n\n")
    cluster_label,finaltime, peak_mem = method_func(adata_raw.copy(),data_name,data_type='Visium',n_clusters=7)
    adata_raw.obs[method_name]=np.array(cluster_label).astype(str)
    result = {
        "method": method_name,
        "exec_time": finaltime,
        "peak_memory": peak_mem
    }
    comp_cost.append(result)
  print("\n\nAll methods have been executed.\n\n")
  return adata_raw, comp_cost
