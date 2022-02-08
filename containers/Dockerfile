FROM docker.io/continuumio/miniconda3

WORKDIR /data/src/

# step 1: build enverioment 
RUN set -xe  \
	&& echo '#!/bin/sh' > /usr/sbin/policy-rc.d  \
	&& echo 'exit 101' >> /usr/sbin/policy-rc.d  \
	&& chmod +x /usr/sbin/policy-rc.d  \
	&& dpkg-divert --local --rename --add /sbin/initctl  \
	&& cp -a /usr/sbin/policy-rc.d /sbin/initctl  \
	&& sed -i 's/^exit.*/exit 0/' /sbin/initctl  \
	&& echo 'force-unsafe-io' > /etc/dpkg/dpkg.cfg.d/docker-apt-speedup  \
	&& echo 'DPkg::Post-Invoke { "rm -f /var/cache/apt/archives/*.deb /var/cache/apt/archives/partial/*.deb /var/cache/apt/*.bin || true"; };' > /etc/apt/apt.conf.d/docker-clean  \
	&& echo 'APT::Update::Post-Invoke { "rm -f /var/cache/apt/archives/*.deb /var/cache/apt/archives/partial/*.deb /var/cache/apt/*.bin || true"; };' >> /etc/apt/apt.conf.d/docker-clean  \
	&& echo 'Dir::Cache::pkgcache ""; Dir::Cache::srcpkgcache "";' >> /etc/apt/apt.conf.d/docker-clean  \
	&& echo 'Acquire::Languages "none";' > /etc/apt/apt.conf.d/docker-no-languages  \
	&& echo 'Acquire::GzipIndexes "true"; Acquire::CompressionTypes::Order:: "gz";' > /etc/apt/apt.conf.d/docker-gzip-indexes  \
	&& echo 'Apt::AutoRemove::SuggestsImportant "false";' > /etc/apt/apt.conf.d/docker-autoremove-suggests
RUN [ -z "$(apt-get indextargets)" ]
RUN mkdir -p /run/systemd  \
	&& echo 'docker' > /run/systemd/container
CMD ["/bin/bash"]

ARG SRC_FOLDER=/data/src
RUN apt-get update -y  \
 && apt install -y gromacs \ 
 && rm -rf /var/lib/apt/lists/*

ENV PATH /opt/conda/bin:$PATH
RUN pip config set global.index-url https://mirrors.aliyun.com/pypi/simple/
RUN conda create -n amber21 -c conda-forge ambertools=21 \
&& rm -rf /opt/conda/pkgs/*

RUN conda init bash \ 
&& echo "conda activate amber21" >> ~/.bashrc

RUN conda install mdanalysis==2.0.0 gmx_MMPBSA pdb2pqr -c conda-forge \
&& rm -rf /opt/conda/pkgs/* 
#RUN conda init bash;conda activate amber21
