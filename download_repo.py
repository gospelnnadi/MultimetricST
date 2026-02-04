

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
