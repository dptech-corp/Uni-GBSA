FROM continuumio/miniconda3

# 1. create an environment
SHELL ["/bin/bash", "-c"]
RUN conda create -n gbsa -c conda-forge acpype openmpi mpi4py gromacs "gmx_MMPBSA>=1.5.6" \
&&  echo 'conda activate gbsa' >> ~/.bashrc \
&&  rm -rf /opt/conda/pkgs/* 

# 2. install unigbsa
RUN source ~/.bashrc \ 
&&  pip install unigbsa lickit \
&&  rm -rf ~/.cache/*

