# Use Ubuntu as base image
FROM ubuntu:22.04

# Set noninteractive for apt
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    wget unzip git curl bzip2 ca-certificates sudo \
    libssl-dev libxml2-dev libcurl4-openssl-dev libopenblas-dev \
    build-essential \
    gfortran \
    libreadline-dev \
    libx11-dev \
    libxt-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libpcre2-dev \
    libjpeg-dev \
    libpng-dev \
    libcairo2-dev \
    libtiff5-dev \
    libblas-dev \
    liblapack-dev \
    && rm -rf /var/lib/apt/lists/*


# Install Miniconda
# Set HOME explicitly (during build this is /root)
ENV HOME=/root

# Install Miniconda like on your local machine
RUN mkdir -p $HOME/miniconda3 && \
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O $HOME/miniconda3/miniconda.sh && \
    bash $HOME/miniconda3/miniconda.sh -b -u -p $HOME/miniconda3 && \
    rm $HOME/miniconda3/miniconda.sh

# Make conda available in PATH
ENV PATH=$HOME/miniconda3/bin:$PATH

RUN conda --version 
# Initialize conda (mainly for interactive shells)
RUN conda init --all

RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main
RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r
RUN conda config --add channels conda-forge


#Create and activate conda environment
RUN conda create -n MMST -y 
SHELL ["conda", "run", "-n", "MMST", "/bin/bash", "-c"]

RUN conda install python=3.10 r-base=4.3.1 r-essentials rpy2=3.5.11 somoclu=1.7.5 -y

ENV R_HOME=/root/miniconda3/envs/MMST/lib/R
ENV LD_LIBRARY_PATH=$R_HOME/lib:$LD_LIBRARY_PATH

# Use bash login shell so conda works
SHELL ["/bin/bash", "--login", "-c"]

# Auto-activate MMST for every shell and process
RUN echo "conda activate MMST" >> ~/.bashrc

# Force all processes to use MMST python first
ENV PATH=/root/miniconda3/envs/MMST/bin:$PATH
ENV CONDA_DEFAULT_ENV=MMST
ENV CONDA_PREFIX=/root/miniconda3/envs/MMST



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

#Install Python requirements
RUN pip install -r requirements.txt

# Install mclust in R
RUN Rscript -e 'install.packages("remotes", repos="https://cran.r-project.org")' && \
    Rscript -e 'remotes::install_version("mclust", version = "6.0.1", repos="https://cran.r-project.org")'

# Expose a working directory for data
VOLUME /workspace

# Clean conda caches
RUN conda clean -a -y && \
    rm -rf /root/miniconda3/pkgs

# Remove pip cache
RUN rm -rf /root/.cache/pip

