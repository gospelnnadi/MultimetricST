
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
    "HERGAST",
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
from runHERGAST import run as runHERGAST 
from runSEDR import run as runSEDR 
from runSpaGCN import run as runSpaGCN
from runSTAGATE import run as runSTAGATE 
from runScanIT import run as runScanIT
from runCCST import run as runCCST
from runConST import run as runConST 
from runSpaceFlow import run as runSpaceFlow 

import pandas as pd
import numpy as np
"""  """
methods = [ 
    ("DeepST", runDeepST) ,
    ("CCST", runCCST),
    ("conST", runConST),
    ("GIST", runGIST), 
    ("GraphST", runGraphST), 
    ("HERGAST", runHERGAST),  
    ("SCAN-IT", runScanIT),  
    ("SEDR", runSEDR),  
    ("SpaceFlow", runSpaceFlow) ,
    ("SpaGCN", runSpaGCN), 
    ("STAGATE", runSTAGATE)
]


import multiprocessing as mp
import numpy as np
import traceback


def _run_single_method(queue, method_name, method_func,
                       adata_raw, data_name, data_type, n_clusters):
    """
    Worker executed in a subprocess.
    """
    try:
        cluster_label, finaltime, peak_mem, allocated, cached = method_func(
            adata_raw,
            data_name,
            data_type=data_type,
            n_clusters=n_clusters
        )

        queue.put({
            "status": "ok",
            "method": method_name,
            "cluster_label": np.array(cluster_label).astype(str),
            "exec_time": finaltime,
            "peak_memory": peak_mem,
            "GPU_exec_time":cached,
            "GPU_peak_memory":  allocated
        })

    except Exception as e:
        queue.put({
            "status": "error",
            "method": method_name,
            "error": str(e),
            "traceback": traceback.format_exc()
        })


def run_clustering_pipeline(adata_raw,
                            data_name,
                            subset_methods=None,
                            data_type='Visium',
                            n_clusters=7,
                            decimal=4):
    
    comp_cost = []

    if subset_methods is None or subset_methods == ["all"]:
        methods_to_run = methods
    else:
        methods_to_run = [(n, f) for n, f in methods if n in subset_methods]

    for method_name, method_func in methods_to_run:

        print(f"\n\nRunning {method_name}...\n\n")

        queue = mp.Queue()

        p = mp.Process(
            target=_run_single_method,
            args=(queue, method_name, method_func,
                  adata_raw.copy(), data_name, data_type, n_clusters)
        )

        p.start()

        try:
            result = queue.get()
        except Exception as e:
            print(f"⚠️ WARNING: {method_name} returned nothing")
            continue

        """  if queue.empty():
            print(f"⚠️ WARNING: {method_name} returned nothing")
            continue 
        if p.exitcode != 0:
            print(f"⚠️ {method_name} crashed (exitcode={p.exitcode})")
            continue"""
        # Skip if child returned None
        if result is None:
            print(f"⚠️ WARNING: {method_name} returned nothing")
            continue

        # Wait for process to finish
        p.join()

        
        if result["status"] == "ok":
            try:
                adata_raw.obs[method_name] = result["cluster_label"]

                comp_cost.append({
                    "method": method_name,
                    "exec_time": result["exec_time"],
                    "peak_memory": result["peak_memory"],
                    "GPU_exec_time": result["GPU_exec_time"],
                    "GPU_peak_memory": result["GPU_peak_memory"]
                })
                
            except Exception as e:
                print(f"⚠️ WARNING: Failed saving results for {method_name}")
                print(e)

        else:
            print(f"⚠️ WARNING: {method_name} failed")
            print(result["error"])
            print(result["traceback"])

        # Force memory cleanup
        del p
        del queue

    print(f"\n\n{subset_methods} methods have been executed.\n\n")

    return adata_raw, comp_cost
