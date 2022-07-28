FROM ubuntu:jammy-20220531
RUN apt update
RUN DEBIAN_FRONTEND=noninteractive apt install -y cmake git g++ doxygen graphviz bash rsync curl mpich libopenblas-dev liblapack-dev liblapacke-dev libeigen3-dev libhdf5-dev libhdf5-mpich-dev clang ccache wget ninja-build wget libarmadillo-dev
RUN export VERSION=5.8.1 && wget https://github.com/GlobalArrays/ga/releases/download/v${VERSION}/ga-${VERSION}.tar.gz && tar xzf ga-${VERSION}.tar.gz && cd ga-${VERSION} && ./configure --disable-f77 --with-mpi3 --without-blas --without-lapack && make install && cd .. && rm -rf ga-${VERSION}*
RUN git clone https://github.com/pjknowles/libtensor && cd libtensor && git checkout gmb && mkdir cmake-build && cd cmake-build && cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=ON -DLIBTENSOR_STANDALONE=ON -DUSE_MKL=0 -G Ninja .. && ninja && ninja install && cd ../.. && rm -rf libtensor
