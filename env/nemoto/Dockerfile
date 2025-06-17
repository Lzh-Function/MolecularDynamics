FROM nvidia/cuda:12.2.2-cudnn8-runtime-ubuntu22.04

WORKDIR /workspace

RUN apt update && \
    apt install -y tzdata

ENV TZ=Asia/Tokyo

RUN apt update \
    && apt install -y \
    curl \
    vim \
    less \
    git \
    nodejs \
    npm \
    wget \
    build-essential \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

RUN apt update \
    && apt install -y locales \
    && locale-gen ja_JP.UTF-8 \
    && echo "export LANG=ja_JP.UTF-8" >> ~/.bashrc

RUN apt update && \
    cd /opt && \
    wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    /bin/bash ./Miniconda3-latest-Linux-x86_64.sh -b -p /opt/miniconda3 && \
    rm ./Miniconda3-latest-Linux-x86_64.sh && \
    export PATH=/opt/miniconda3/bin:$PATH && \
    echo ". /opt/miniconda3/etc/profile.d/conda.sh" >> ~/.bashrc && \
    ln -s ./miniconda3/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo "conda activate" >> ~/.bashrc && \
    . /opt/miniconda3/etc/profile.d/conda.sh && \
    conda activate && \
    conda update conda && \
    conda create -y -n mdenv python=3.11 && \
    conda activate mdenv && \
    conda install -y -c conda-forge rdkit jupyter ipykernel ipywidgets nodejs ncurses ipympl notebook jupyter_contrib_nbextensions && \
    pip install jupyterlab-widgets && \
    conda install -y -c conda-forge openmm cudatoolkit=11.8 && \
    conda install -y -c conda-forge openmmforcefields openff-toolkit pdbfixer parmed mdtraj mdanalysis nglview pdb2pqr && \
    git clone https://github.com/openbabel/openbabel.git && \
    apt install -y libxml2-dev libboost-all-dev libbz2-dev libomp-dev zlib1g-dev libeigen3-dev libcairo2-dev cmake swig && \
    cd openbabel && \
    mkdir build && cd build && \
    cmake -DENABLE_OPENMP=ON -DBUILD_GUI=OFF -DPYTHON_EXECUTABLE=/opt/miniconda3/envs/mdenv/bin/python3 -DPYTHON_BINDINGS=ON -DRUN_SWIG=ON .. && \
    make -j 33 && \
    make install && \
    conda config --add channels conda-forge && \
    pip install git+https://github.com/chemosim-lab/ProLIF.git && \
    rm -rf /var/lib/apt/lists/* 

ENV LD_LIBRARY_PATH=/opt/openbabel/build/lib:${LD_LIBRARY_PATH} \
    PATH=/opt/openbabel/build/bin:${PATH}