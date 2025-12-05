#!/bin/bash

module load cuda11.4/toolkit/11.4.1  
module load openmpi/cuda/64/3.1.1 

. "/home/accounts/personale/nndgpl46/miniconda3/etc/profile.d/conda.sh"

conda activate bertwalk2 #mmst 

cd /home/accounts/personale/nndgpl46/MultimetricST/

python MultimetricST.py  --data_path $1  --ground_truth $2 --mode $3 --is_h5ad $4 --method_cluster_label $5 $6 --result_savepath $7  


