# Use Ubuntu as base image
FROM ubuntu:22.04

# Set noninteractive for apt
ENV DEBIAN_FRONTEND=noninteractive


# Install Miniconda
ENV CONDA_DIR=/opt/conda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-py310_23.11.1-0-Linux-x86_64.sh -O /tmp/miniconda.sh && \
    bash /tmp/miniconda.sh -b -p $CONDA_DIR && \
    rm /tmp/miniconda.sh
ENV PATH=$CONDA_DIR/bin:$PATH

# Create and activate conda environment
RUN conda create -n MMST -y python=3.10 r-base=4.3.1 somoclu=1.7.5
SHELL ["conda", "run", "-n", "MMST", "/bin/bash", "-c"]

# Install PyTorch and PyTorch Geometric
RUN pip install torch==2.1.0 \
    torchaudio==2.1.0 \
    torchvision==0.16.0

RUN pip install torch_sparse -f https://data.pyg.org/whl/torch-2.1.0+cu121.html \
    torch_scatter -f https://data.pyg.org/whl/torch-2.1.0+cu121.html \
    torch_cluster -f https://data.pyg.org/whl/torch-2.1.0+cu121.html \
    torch_spline_conv -f https://data.pyg.org/whl/torch-2.1.0+cu121.html

RUN pip install torch-geometric==2.7.0

# Copy MultimetricST code into the container
WORKDIR /workspace
COPY . /workspace

# Install Python requirements
RUN pip install -r requirements.txt

# Install mclust in R
RUN Rscript -e 'install.packages("remotes", repos="https://cran.r-project.org")' && \
    Rscript -e 'remotes::install_version("mclust", version = "6.0.1", repos="https://cran.r-project.org")'

# Expose a working directory for data
VOLUME /workspace

# Default command to run MultimetricST
# Can be overridden with `docker run ...`
#CMD ["conda", "run", "-n", "MMST", "python", "MultimetricST.py", "--mode", "1", "--data_path", "Data/DLPFC/151673", "--ground_truth", "Data/DLPFC/151673/metadata.tsv", "--ground_truth_col_name", "layer_guess", "--data_name", "151673", "--is_h5ad", "0", "--data_type", "Visium", "--n_clusters", "7"]
