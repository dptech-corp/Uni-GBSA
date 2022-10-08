FROM continuumio/miniconda3

# 1. create a enverioment
SHELL ["/bin/bash", "-c"]
RUN conda create -n gbsa -c conda-forge acpype gromacs openmpi mpi4py \
&&  echo -e 'conda activate gbsa' >> ~/.bashrc \
&&  echo -e 'alias ls="ls --color" \n' >> ~/.bashrc \
&&  rm -rf /opt/conda/pkgs/* 


# 2. install Uni-GBSA
RUN source ~/.bashrc \ 
&&  pip install unigbsa gmx_MMPBSA>=1.5.6 lickit \
&&  rm -rf ~/.cache/*

