
---

# READMe for Containerized (Docker) Environment

## Setup

We provide a **Docker-based environment** to ensure full reproducibility and compatibility with the evaluated methods.

install MultimetricST package.
````
git clone https://github.com/InfOmics/MultimetricST.git
cd ~/MultimetricST
````
To use the mmst container install [Docker](https://www.docker.com/) (version>=28.4.0 is recommend) and [Docker Compose] (https://docs.docker.com/compose/).

Download packages of the spatial transcriptomics spatial domain identification methods to be evaluated described in the paper.
````
python download_repo.py
````

If python not installed. 
````
sudo apt install python3-pip
````
````
python3 download_repo.py 
````
Repositories will be stored in:

MultimetricST/Spatial_Clustering_Methods/

## Start the MultimetricST container
Get the mmst container with docker compose.
`````
docker-compose up -d
`````
Verify the container image.
````` 
docker images
`````

### Data Availability ###
The spatial transcriptomics datasets are available at:  https://zenodo.org/records/18482658

Download the DLPFC 151673 data:

        wget https://zenodo.org/records/18482658/files/Data.zip
        unzip Data.zip
        rm Data.zip

Download the Axolotl dataset:

        wget https://zenodo.org/records/18482658/files/Stereo.zip
        unzip Stereo.zip -d Data/
        rm Stereo.zip

if unzip not installed.
````
sudo apt-get install zip unzip
````

### Access the container environment.
`````
docker-compose exec mmst bash
`````
For detailed command-line usage and dataset-specific examples, see [View Usage Documentation](USAGE.md).

To exit the container environment. 
`````
Ctrl+D
`````
To stop the container.
`````
docker-compose down
`````


