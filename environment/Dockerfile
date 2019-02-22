FROM ubuntu:bionic
MAINTAINER Björn Dahlgren <bjodah@DELETEMEgmail.com>
ENV LANG=en_US.UTF-8 LANGUAGE=en_US:en LC_ALL=en_US.UTF-8 DEBIAN_FRONTEND=noninteractive DEBCONF_NONINTERACTIVE_SEEN=true CMAKE_PREFIX_PATH=/usr/lib/llvm-7:/opt/syme-af64b72 CC=gcc-8 CXX=g++-8 FC=gfortran-8
RUN apt-get update && \
    apt-get --quiet --assume-yes install locales && \
    locale-gen en_US.UTF-8 && \
    echo "path-exclude /usr/share/doc/*" >/etc/dpkg/dpkg.cfg.d/01_nodoc && \
    echo "path-include /usr/share/doc/*/copyright" >>/etc/dpkg/dpkg.cfg.d/01_nodoc && \
    apt-get upgrade -y && \
    apt-get --assume-yes --no-install-recommends install curl ssh gnupg2 time sudo git zlib1g-dev \
    gcc-8 g++-8 gfortran-8 \
    make cmake binutils-dev bzip2 libgmp-dev libgsl-dev \
    liblapack-dev python3-dev python3-pip python3-setuptools python3-wheel && \
    sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-8 800 --slave /usr/bin/g++ g++ /usr/bin/g++-8 && \
    curl -Ls https://apt.llvm.org/llvm-snapshot.gpg.key|sudo apt-key add - && \
    echo "deb http://apt.llvm.org/bionic/ llvm-toolchain-bionic-7 main" | tee -a /etc/apt/sources.list && \
    apt-get update && apt-get --quiet --assume-yes --no-install-recommends install llvm-7-dev && \
    python3 -m pip install --upgrade pip setuptools wheel && \
    python3 -m pip install http://hera.physchem.kth.se/~repo/616a6ebc248f81b461c0562aa942de37/nbconvert-ad27a93.tar.gz && \
    python3 -m pip install future appdirs pybind11 cython pytest numpy scipy statsmodels "sympy>=1.1.1,!=1.2" matplotlib jupyter ipython ipykernel \
            scikit-optimize quantities pyparsing ipywidgets pulp wurlitzer argh mako && \
    python3 -m ipykernel install && \
    python3 -m jupyter nbextension enable --py widgetsnbextension && \
    echo "backend: Agg" > /etc/matplotlibrc && \
    mkdir -p /root/.config/matplotlib/ && cp /etc/matplotlibrc /root/.config/matplotlib/ && \
    python3 -c "import matplotlib.pyplot as plt" && \
    apt-get clean && \
    cd /tmp && curl -Ls https://dl.bintray.com/boostorg/release/1.69.0/source/boost_1_69_0.tar.bz2 | tar xj && \
    cd boost_1_69_0 && echo $CC && ./bootstrap.sh --with-toolset=gcc --prefix=/usr || cat bootstrap.log && ./b2 install >b2_build.log && tail b2_build.log && \
    ln -s $(which gcc-8) /usr/bin/cc && \
    ln -sf $(which gcc-8) /usr/bin/x86_64-linux-gnu-gcc && \
    ln -sf $(which g++-8) /usr/bin/x86_64-linux-gnu-g++ && \
    curl -Ls http://computation.llnl.gov/projects/sundials/download/sundials-3.2.1.tar.gz | tar xz && \
    mkdir build/ && cd build/ && \
    cmake -DCMAKE_INSTALL_PREFIX:PATH=/usr -DCMAKE_BUILD_TYPE:STRING="Release" -DBUILD_SHARED_LIBS:BOOL=ON -DBUILD_STATIC_LIBS:BOOL=OFF \
    -DEXAMPLES_ENABLE_C:BOOL=OFF -DEXAMPLES_INSTALL:BOOL=OFF -DOPENMP_ENABLE:BOOL=OFF -DLAPACK_ENABLE:BOOL=ON -DSUNDIALS_INDEX_SIZE=32 \
    ../sundials*/ && make install && cd - && rm -r build/ sundials* && \
    mkdir /tmp/syme && cd /tmp/syme && curl -Ls https://github.com/symengine/symengine/archive/af64b726ea939f206380107b3855dd7ad277fe71.tar.gz | tar xz && \
    cd /tmp/syme && mkdir build && cd build && cmake -DCMAKE_INSTALL_PREFIX:PATH=/opt/syme-af64b72 -DCMAKE_BUILD_TYPE:STRING="Release" -DWITH_BFD:BOOL=ON -DBUILD_SHARED_LIBS:BOOL=ON \
    -DWITH_LLVM:BOOL=ON -DBUILD_TESTS:BOOL=OFF -DBUILD_BENCHMARKS:BOOL=OFF ../symengine-af64b726ea939f206380107b3855dd7ad277fe71 && \
    make install && cd .. && rm -r build symengine-*/ && \
    curl -Ls https://github.com/symengine/symengine.py/archive/ab5841e26f87aaa3e540b14c985b2133d97de07f.tar.gz | tar xz && \
    cd symengine.py-ab5841e26f87aaa3e540b14c985b2133d97de07f && python3 setup.py install && cd .. && rm -r symengine.py-* && \
    python3 -c 'import symengine; assert symengine.test()' && \
    python3 -m pip install --upgrade "pyodeint>=0.9.8,<0.10.0" "pygslodeiv2>=0.8.4,<0.9.0" "pycvodes>=0.10.13,<0.11.0" "pycodeexport>=0.1.2,<0.2.0" "pycompilation>=0.4.3,<0.5.0" notebook "matplotlib>=2.2.3,<3.0.0" "pykinsol>=0.1.3" && \
    python3 -m pytest --pyargs pyodeint && \
    python3 -m pytest --pyargs pygslodeiv2 && \
    python3 -m pytest --pyargs pycvodes && \
    python3 -m pytest --pyargs pykinsol && \
    python3 -m pip install --upgrade "pyneqsys>=0.5.6,<0.6.0" "pyodesys>=0.11.17,<0.12.0" "sympy==1.3" "block_diag_ilu>=0.3.9,<0.4.0" "finitediff>=0.6.2" && \
    python3 -m pytest --pyargs pyneqsys && \
    python3 -m pip install "chempy>=0.6.16,<0.7.0" && \
    python3 -m pytest --pyargs chempy && \
    python3 -m pip install "chemreac>=0.7.15,<0.8.0" && \
    python3 -m pytest --pyargs chemreac && \
    curl -Ls http://hera.physchem.kth.se/~kinetgine/releases/pykinetgine-0.3.5.tar.gz | tar xz -C /opt && \
    python3 -m pip install http://hera.physchem.kth.se/~aq_radiolysis/releases/aq_radiolysis-0.1.4.tar.gz && \
    python3 -m pytest --pyargs aq_radiolysis && \
    python3 -m pip install http://hera.physchem.kth.se/~chemreac-derivations/releases/advediffreac-0.1.3.tar.gz && \
    python3 -m pytest --pyargs advediffreac && \
    python3 -m pip install http://hera.physchem.kth.se/~batemaneq/releases/batemaneq-0.2.0.tar.gz && \
    python3 -m pytest --pyargs batemaneq && \    
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* ~/.cache/pip
