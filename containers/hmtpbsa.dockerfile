FROM continuumio/miniconda3
# install
RUN apt-get update -y && \
	apt-get install gromacs -y && \ 
	apt-get install vim -y && \
        rm -rf /var/lib/apt/lists/*
# python env
RUN conda create -n amber21 python=3.9 \
                            ambertools=21 \
                            pip=22.0.3 \
                            mdanalysis=2.0.0 \
                            gmx_MMPBSA=1.5.2 \
                            pdb2pqr=3.4.1 \
                            pytest=6.2.5 \
                            pytest-cov=3.0.0 \
                            pylint=2.12.2 \
                            acpype=2021.02 \
                            -c conda-forge -y && \
                            # clean up
                            rm -rf /opt/conda/pkgs/*
                            
# add env
RUN echo "source activate amber21" >> ~/.bashrc
