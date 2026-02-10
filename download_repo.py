

import sys, os
ROOT = os.getcwd()
print("The current working dir is, ", ROOT)


import subprocess
import os

def download_repo():
    target_dir = f"{ROOT}/Spatial_Clustering_Methods"
    os.makedirs(target_dir, exist_ok=True)  # ensure target directory exists

    repos_file = os.path.join(target_dir, "repos_git.txt")

    # Read repo URLs from file (one per line)
    with open(repos_file, "r") as f:
        repos = [line.strip() for line in f if line.strip()]

    for repo in repos:
        repo_name = os.path.splitext(os.path.basename(repo))[0]
        repo_path = os.path.join(target_dir, repo_name)

        if os.path.exists(repo_path):
            print(f"✅ Skipping {repo_name} (already exists)")
            continue

        print(f"⬇️ Cloning {repo} into {repo_path} ...")
        subprocess.run(["git", "clone", "--depth", "1", repo, repo_path], check=True)
    
    print("✨ All repositories processed!")

download_repo()

def fix_known_incompactibilities():
    from pathlib import Path

    # --- SPACEFLOW FIX ---
    spaceflow_file = Path("Spatial_Clustering_Methods/SpaceFlow/SpaceFlow/SpaceFlow.py")
    if spaceflow_file.exists():
      
        text = spaceflow_file.read_text()

        old_line = "sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, flavor='cell_ranger', subset=True)"
        new_line = "sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, subset=True)"

        if old_line in text:
            text = text.replace(old_line, new_line)
            spaceflow_file.write_text(text)
            print("✅  SpaceFlow patch applied")
        else:
            print("SpaceFlow line not found — already patched")


    # --- SEDR FIX ---
    sedr_file = Path("Spatial_Clustering_Methods/SEDR/SEDR/clustering_func.py")
    if sedr_file.exists():
        lines = sedr_file.read_text().splitlines()

        with sedr_file.open("w") as f:
            for line in lines:
                # Comment out any R home forcing lines
                if "R_HOME" in line and "#" not in line:
                    f.write("\t#" + line + "\n")
                    print("✅ SEDR R-environment patch applied")
                else:
                    f.write(line + "\n")
   

fix_known_incompactibilities()