FROM --platform=linux/amd64 ubuntu:26.04

RUN apt-get update && \
      DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
      make \
      meson \
      cmake \
      g++ \
      gfortran \
      git \
      htop \
      nano \
      vim \
      less \
      curl \
      hdf5-tools \
      openssh-client \
      patch \
      pkg-config \
      libboost-dev \
      libboost-mpi-dev \
      libboost-serialization-dev \
      libeigen3-dev \
      libfftw3-dev \
      libgfortran5 \
      libgmp-dev \
      libgsl-dev \
      libmpfr-dev \
      libhdf5-dev \
      libmkl-dev \
      libopenmpi-dev \
      libnfft3-dev \
      libzstd-dev \
      openmpi-bin \
      openmpi-common \
      openmpi-doc \
      python3-dev \
      python3-h5py \
      python3-ipython \
      python3-mako \
      python3-matplotlib \
      python3-mpi4py \
      python3-numpy \
      python3-pandas \
      python3-pip \
      python3-pytest \
      python3-scipy \
      python3-setuptools \
      python3-skimage \
      python3-tk \
      python3-tomli \
      python3-venv \
      python3-yaml \
      && \
    apt-get autoremove --purge -y && \
    apt-get autoclean -y && \
    rm -rf /var/cache/apt/* /var/lib/apt/lists/*

RUN ln -s /usr/bin/python3 /usr/bin/python

ARG NB_USER=triqs
ARG NB_UID=1000
RUN useradd -u ${NB_UID} -m ${NB_USER} -o

USER ${NB_USER}
RUN python3 -m venv --system-site-packages /home/${NB_USER}/.venv
ENV VIRTUAL_ENV=/home/${NB_USER}/.venv \
    PATH=/home/${NB_USER}/.venv/bin:$PATH

ARG BRANCH=unstable
ARG NCORES=10
ARG ARCH=native
ENV SRC=/tmp/src BUILD=/tmp/build INSTALL=/usr \
    CPLUS_INCLUDE_PATH=/usr/include/x86_64-linux-gnu/openmpi:/usr/include/hdf5/serial:/usr/include/mkl:$CPLUS_INCLUDE_PATH \
    C_INCLUDE_PATH=/usr/include/python3.14 \
    CC=gcc CXX=g++ CXXFLAGS="-march=${ARCH}"

# -- TRIQS ecosystem (cmake) --------------------------------------------------
USER root
RUN set -ex ; \
  for pkg in triqs cthyb ctseg ctint ; do \
    mkdir $BUILD ; cd $BUILD ; \
    git clone https://github.com/TRIQS/$pkg --branch DLR2D --depth 1 $SRC ; \
    cmake $SRC -DCMAKE_INSTALL_PREFIX=$INSTALL ; \
    make -j$NCORES ; \
    make install ; \
    rm -rf $SRC $BUILD ; \
  done

RUN set -ex ; \
  for pkg in tprf ; do \
    mkdir $BUILD ; cd $BUILD ; \
    git clone https://github.com/TRIQS/$pkg --branch $BRANCH --depth 1 $SRC ; \
    cmake $SRC -DCMAKE_INSTALL_PREFIX=$INSTALL ; \
    make -j$NCORES ; \
    make install ; \
    rm -rf $SRC $BUILD ; \
  done

# -- pyed (pure Python ED) ----------------------------------------------------
USER ${NB_USER}
RUN pip install git+https://github.com/HugoStrand/pyed.git

# -- libcommute + pomerol + pomerol2triqs (cmake) -----------------------------
USER root
RUN set -ex ; \
  mkdir $BUILD ; cd $BUILD ; \
  git clone https://github.com/krivenko/libcommute --depth 1 $SRC ; \
  cmake $SRC -DCMAKE_INSTALL_PREFIX=$INSTALL -DTESTS=OFF -DEXAMPLES=OFF ; \
  make install ; \
  rm -rf $SRC $BUILD

RUN set -ex ; \
  mkdir $BUILD ; cd $BUILD ; \
  git clone https://github.com/aeantipov/pomerol --depth 1 $SRC ; \
  cmake $SRC -DCMAKE_INSTALL_PREFIX=$INSTALL -DTesting=OFF ; \
  make -j$NCORES ; \
  make install ; \
  rm -rf $SRC $BUILD

RUN set -ex ; \
  mkdir $BUILD ; cd $BUILD ; \
  git clone https://github.com/krivenko/pomerol2triqs --depth 1 $SRC ; \
  cmake $SRC -DCMAKE_INSTALL_PREFIX=$INSTALL ; \
  make -j$NCORES ; \
  make install ; \
  rm -rf $SRC $BUILD

# -- SciFortran + EDIpack + edipack2py + edipack2triqs -------------------------
RUN set -ex ; \
  mkdir $BUILD ; cd $BUILD ; \
  git clone https://github.com/SciFortran/SciFortran --depth 1 $SRC ; \
  cmake $SRC -DCMAKE_INSTALL_PREFIX=$INSTALL ; \
  make -j$NCORES ; \
  make install ; \
  rm -rf $SRC $BUILD

ENV PKG_CONFIG_PATH=$INSTALL/etc:$PKG_CONFIG_PATH
RUN set -ex ; \
  mkdir $BUILD ; cd $BUILD ; \
  git clone https://github.com/EDIpack/EDIpack --depth 1 $SRC ; \
  cmake $SRC -DCMAKE_INSTALL_PREFIX=$INSTALL ; \
  make -j$NCORES ; \
  make install ; \
  rm -rf $SRC $BUILD

USER ${NB_USER}
RUN pip install edipack2py
RUN pip install git+https://github.com/krivenko/edipack2triqs.git

# -- w2dynamics_interface (cmake, fetches w2dynamics via FetchContent) ----------
USER ${NB_USER}
RUN pip install configobj
USER root
RUN set -ex ; \
  mkdir $BUILD ; cd $BUILD ; \
  git clone https://github.com/TRIQS/w2dynamics_interface --branch $BRANCH --depth 1 $SRC ; \
  cmake $SRC -DCMAKE_INSTALL_PREFIX=$INSTALL -DBuild_Tests=OFF ; \
  sed -i '/ADD_SUBDIRECTORY(testsuite/d; /enable_testing/d' $BUILD/_deps/w2dynamics-src/CMakeLists.txt ; \
  cmake $SRC ; \
  make -j$NCORES ; \
  make install ; \
  rm -rf $SRC $BUILD

# -- nrgljubljana_interface (cmake, fetches nrgljubljana via FetchContent) -----
RUN set -ex ; \
  mkdir $BUILD ; cd $BUILD ; \
  git clone https://github.com/TRIQS/nrgljubljana_interface --branch $BRANCH --depth 1 $SRC ; \
  cmake $SRC -DCMAKE_INSTALL_PREFIX=$INSTALL ; \
  make -j$NCORES ; \
  make install ; \
  rm -rf $SRC $BUILD

# -- forktps (cmake, private repo -- skipped) ----------------------------------
# RUN set -ex ; \
#   mkdir $BUILD ; cd $BUILD ; \
#   git clone https://github.com/TRIQS/forktps --branch $BRANCH --depth 1 $SRC ; \
#   cmake $SRC -DCMAKE_INSTALL_PREFIX=$INSTALL ; \
#   make -j$NCORES ; \
#   make install ; \
#   rm -rf $SRC $BUILD

# -- ALPSCore + ALPS/CT-HYB (cmake) --------------------------------------------
RUN set -ex ; \
  mkdir $BUILD ; cd $BUILD ; \
  git clone https://github.com/ALPSCore/ALPSCore --depth 1 $SRC ; \
  cmake $SRC -DCMAKE_INSTALL_PREFIX=$INSTALL -DCMAKE_POLICY_VERSION_MINIMUM=3.5 -DENABLE_MPI=ON -DTesting=OFF ; \
  make -j$NCORES ; \
  make install ; \
  rm -rf $SRC $BUILD

RUN set -ex ; \
  mkdir $BUILD ; cd $BUILD ; \
  git clone https://github.com/ALPSCore/CT-HYB --depth 1 $SRC ; \
  cmake $SRC -DCMAKE_INSTALL_PREFIX=$INSTALL -DCMAKE_POLICY_VERSION_MINIMUM=3.5 -DTesting=OFF ; \
  make -j$NCORES ; \
  make install ; \
  rm -rf $SRC $BUILD

# -- DCore + dcorelib (pip) ----------------------------------------------------
USER ${NB_USER}
RUN pip install 'sparse_ir<2' cvxpy toml sympy
RUN pip install --no-deps git+https://github.com/shinaoka/dcorelib.git && \
    python -c "import dcorelib, pathlib; \
      f = pathlib.Path(dcorelib.__file__); \
      f.write_text(f.read_text() + '\nfrom importlib.metadata import version as _v\n__version__ = _v(\"dcorelib\")\n')"
RUN pip install --no-deps git+https://github.com/issp-center-dev/DCore.git@develop

# -- Runtime setup -------------------------------------------------------------
USER ${NB_USER}
WORKDIR /home/${NB_USER}/benchmarks
ENV OMPI_ALLOW_RUN_AS_ROOT=1 \
    OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 \
    OMPI_MCA_btl_vader_single_copy_mechanism=none \
    MKL_THREADING_LAYER=SEQUENTIAL

CMD ["/bin/bash"]
