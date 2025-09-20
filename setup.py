from setuptools import setup, find_packages

with open("requirements.txt") as f:
    requirements = f.read().splitlines()


setup(
    name='MultimetricST',
    version='0.1.0',
    author='Gospel Ozioma Nnadi',
    description='A Multi-Perspective Evaluation Framework of Spatial Transcriptomics Clustering Methods',
    author_email ='gospelozioma.nnadi@univr.it',
    license = 'MIT',
    packages =find_packages(), 
    install_requires=requirements,
    python_requires='>=3.10',
    zip_safe = False,
    include_package_data = True,
    url='https://github.com/InfOmics/MultimetricST.git',
)