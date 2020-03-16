FROM ubuntu:xenial

# metadata
LABEL base.image="ubuntu:xenial"
LABEL container.version="1"
LABEL software="ARTIC-nCov19"
LABEL software.version="1.0"
LABEL description="Conda environment for ARTIC network nCov19 bioinformatic SOP"
LABEL website="https://artic.network/ncov-2019/ncov2019-bioinformatics-sop.html"
LABEL license=""
LABEL maintainer="Kelsey Florek"
LABEL maintainer.email="kelsey.florek@slh.wisc.edu"

# install needed software for conda install
RUN apt-get update && apt-get install -y \
  wget \
  git \
  build-essential

# get miniconda and the artic-ncov19 repo
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh &&\
 bash ./Miniconda3-latest-Linux-x86_64.sh -p /miniconda -b  &&\
  rm Miniconda3-latest-Linux-x86_64.sh &&\
  git clone --recursive https://github.com/artic-network/artic-ncov2019.git

# set the environment
ENV PATH="/miniconda/bin:$PATH"

# create the conda environment
RUN conda env create -f /artic-ncov2019/environment.yml
RUN conda init
WORKDIR /data
