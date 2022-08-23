FROM continuumio/miniconda3

# 1. create a enverioment
SHELL ["/bin/bash", "-c"]
RUN conda create -n amber21 -c conda-forge ambertools=21 acpype=2021.02 gromacs openmpi mpi4py \
&&  echo 'conda activate amber21' >> ~/.bashrc \
&&  rm -rf /opt/conda/pkgs/* 


# 2. install requirments
RUN source ~/.bashrc \ 
&&  pip install gmx_MMPBSA \
&&  rm -rf ~/.cache/*

