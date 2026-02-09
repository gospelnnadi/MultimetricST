
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
"""   """
methods = [ 
    ("DeepST", runDeepST) ,
    ("CCST", runCCST),
    ("conST", runConST),
    ("GIST", runGIST), 
    ("GraphST", runGraphST), 
    ("SCAN-IT", runScanIT), 
    ("SEDR", runSEDR),
    ("SpaceFlow", runSpaceFlow) ,
    ("SpaGCN", runSpaGCN), 
    ("STAGATE", runSTAGATE)
]
comp_cost = []
def run_clustering_pipeline(adata_raw,data_name,subset_methods=None,data_type='Visium',n_clusters=7,decimal=4):
  comp_cost = []
  print("show subset_methods 2",  subset_methods)

  if subset_methods is None or subset_methods == ["all"]:
        methods_to_run = methods
  else:
        print("show subset_methods3 ", subset_methods)
        methods_to_run = [(n, f) for n, f in methods if n in subset_methods]
        print(f"Selected methods: {[n for n, _ in methods_to_run]}")
  for method_name, method_func in methods_to_run:
    
    print(f"\n\nRunning {method_name}...\n\n\n\n\n")
    cluster_label,finaltime, peak_mem = method_func(adata_raw.copy(),data_name,data_type=data_type,n_clusters=7)
    adata_raw.obs[method_name]=np.array(cluster_label).astype(str)
    result = {
        "method": method_name,
        "exec_time": finaltime,
        "peak_memory": peak_mem
    }
    comp_cost.append(result)
  print(f"\n\n{subset_methods} methods have been executed.\n\n")
  return adata_raw, comp_cost
