FROM nvidia/cuda:12.4.1-cudnn-devel-ubuntu22.04

ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=Asia/Tokyo

ARG USERNAME=function
ARG USER_UID=1000
ARG USER_GID=${USER_UID}

WORKDIR /workspace

RUN apt-get update && apt-get install -y --no-install-recommends \
    curl vim less git tzdata nodejs npm \
    htop wget sudo pv bash-completion openssh-client \
    make build-essential gfortran flex bison bc openmpi-bin libssl-dev zlib1g-dev \
    libbz2-dev libreadline-dev libsqlite3-dev \
    libxml2-dev libboost-all-dev libomp-dev libfftw3-dev libnetcdf-dev libblas-dev liblapack-dev \
    zlib1g-dev libeigen3-dev libcairo2-dev swig lsof \
    libncursesw5-dev xz-utils tk-dev libxmlsec1-dev libffi-dev liblzma-dev \
    libglm-dev libglew-dev openbabel libopenbabel-dev software-properties-common lsb-release && \
    wget -O - https://apt.kitware.com/keys/kitware-archive-latest.asc 2>/dev/null | gpg --dearmor - | tee /etc/apt/trusted.gpg.d/kitware.gpg >/dev/null && \
    apt-add-repository -y "deb https://apt.kitware.com/ubuntu/ $(lsb_release -cs) main" && \
    apt-get update && \
    apt-get install -y cmake && \
    add-apt-repository ppa:openslide/openslide && \
    apt-get install -y openslide-tools python3-openslide && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# add USER
RUN groupadd --gid ${USER_GID} ${USERNAME} && \
    useradd --uid ${USER_UID} --gid ${USER_GID} -m -s /bin/bash ${USERNAME} && \
    echo "${USERNAME} ALL=(ALL) NOPASSWD: ALL" >> /etc/sudoers && \
    chown -R ${USERNAME}:${USERNAME} /workspace && \
    chown -R ${USERNAME}:${USERNAME} /opt

USER ${USERNAME}

# OpenMM Python Environment (conda)
RUN cd /opt && \
    wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    /bin/bash ./Miniconda3-latest-Linux-x86_64.sh -b -p /opt/miniconda3 && \
    rm ./Miniconda3-latest-Linux-x86_64.sh && \
    sudo ln -s ./miniconda3/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/miniconda3/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate vsenv_830" >> ~/.bashrc

ENV PATH=/opt/miniconda3/bin:${PATH}

RUN conda update -n base -c defaults conda -y && \
    # openmm 8.3.0-beta
    conda create -y -n vsenv_830 python=3.11.10 && \
    conda run -n vsenv_830 conda install -y -c conda-forge jupyter ipykernel ipywidgets nodejs ncurses ipympl notebook jupyter_contrib_nbextensions mamba \
    gcc_linux-64 gxx_linux-64 && \
    conda run -n vsenv_830 mamba install -y -c conda-forge numpy pandas matplotlib seaborn scipy tqdm && \
    conda run -n vsenv_830 mamba install -y -c conda-forge/label/openmm_rc -c conda-forge rdkit openmm openmmforcefields openff-toolkit pyscf && \
    conda run -n vsenv_830 mamba install -y -c conda-forge pdbfixer parmed mdtraj==1.10.3 mdanalysis nglview pdb2pqr plip prolif pymol-open-source geometric py3Dmol && \
    conda clean -afy && \
    # openmm 8.2.0
    conda create -y -n vsenv_820 python=3.11.10 && \
    conda run -n vsenv_820 conda install -y -c conda-forge jupyter ipykernel ipywidgets nodejs ncurses ipympl notebook jupyter_contrib_nbextensions mamba \
    gcc_linux-64 gxx_linux-64 && \
    conda run -n vsenv_820 mamba install -y -c conda-forge numpy pandas matplotlib seaborn scipy tqdm && \
    conda run -n vsenv_820 mamba install -y -c conda-forge rdkit openmm==8.2.0 openmmforcefields openff-toolkit pyscf && \
    conda run -n vsenv_820 mamba install -y -c conda-forge pdbfixer parmed mdtraj==1.10.3 mdanalysis nglview pdb2pqr plip prolif pymol-open-source geometric py3Dmol && \
    conda clean -afy

# Amber24(pmemd24)
COPY pmemd24.tar.bz2 /opt/

RUN cd /opt && \
    tar xvfj pmemd24.tar.bz2 && \
    rm pmemd24.tar.bz2 && \
    cd pmemd24_src && \
    ./update_amber --update-to ambertools/24.10,Amber/24.1 && \
    cd build && \
    sed -i s/-DMPI=FALSE\ -DCUDA=FALSE/-DMPI=TRUE\ -DOPENMP=TRUE\ -DCUDA=TRUE/g run_cmake && \
    ./run_cmake && \
    make -j16 install

ENV PMEMDHOME=/opt/pmemd24
ENV PATH="${PMEMDHOME}/bin:${PATH}"

# AmberTools25
COPY ambertools25.tar.bz2 /opt/

RUN cd /opt && \
    tar xvfj ambertools25.tar.bz2 && \
    rm ambertools25.tar.bz2 && \
    cd ambertools25_src && \
    ./update_amber --update && \
    cd build && \
    sed -i s/-DMPI=FALSE\ -DCUDA=FALSE/-DMPI=TRUE\ -DOPENMP=TRUE\ -DCUDA=TRUE/g run_cmake && \
    sed -i s/-DDOWNLOAD_MINICONDA=TRUE/-DDOWNLOAD_MINICONDA=TRUE\ -DCMAKE_POLICY_VERSION_MINIMUM=3.5/g run_cmake && \
    ./run_cmake && \
    make -j16 install 

ENV AMBERHOME=/opt/ambertools25
ENV PATH="${AMBERHOME}/bin:${PATH}"
ENV LD_LIBRARY_PATH="${AMBERHOME}/lib:${LD_LIBRARY_PATH}"
ENV PYTHONPATH=/opt/ambertools25/lib/python3.12/site-packages
ENV PERL5LIB=/opt/ambertools25/lib/perl

# GROMACS
COPY gromacs-2024.2.tar.gz /opt/

RUN cd /opt && \
    tar zxf gromacs-2024.2.tar.gz && \
    rm gromacs-2024.2.tar.gz && \
    cd gromacs-2024.2 && \
    mkdir build && \
    cd build && \
    cmake .. -DGMX_BUILD_OWN_FFTW=ON -DREGRESSIONTEST_DOWNLOAD=ON -DCMAKE_INSTALL_PREFIX=/opt/gromacs-2024.2 \
    -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DGMX_MPI=OFF -DGMX_DOUBLE=OFF -DGMX_OPENMP=ON -DGMX_GPU=CUDA && \
    make -j16 && \
    make check && \
    make install && \
    . /opt/gromacs-2024.2/bin/GMXRC && \
    echo "source /opt/gromacs-2024.2/bin/GMXRC" >> ~/.bashrc

# Rosetta
COPY rosetta_src_3.14_bundle.tar.bz2 /opt/

RUN cd /opt && \
    pv rosetta_src_3.14_bundle.tar.bz2 | tar -jxf - && \
    rm rosetta_src_3.14_bundle.tar.bz2 && \
    mv rosetta.source.release-371 rosetta_3.14 && \
    cd rosetta_3.14/main/source && \
    conda run -n vsenv_830 ./scons.py -j16 mode=release bin && \
    conda run -n vsenv_830 ./scons.py -j16 bin mode=release extras=mpi

CMD ["/bin/bash"]