#FROM ubuntu:groovy-20200812
FROM intel/oneapi-hpckit

RUN apt-get update 
RUN apt-get upgrade -y
#RUN DEBIAN_FRONTEND=noninteractive apt-get install -y cmake git g++ gfortran doxygen graphviz bash rsync curl
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y cmake git doxygen graphviz rsync curl
#RUN DEBIAN_FRONTEND=noninteractive apt-get install -y mpich
#RUN DEBIAN_FRONTEND=noninteractive apt-get install -y libblas-dev liblapack-dev
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y libeigen3-dev
#RUN DEBIAN_FRONTEND=noninteractive apt-get install -y libhdf5-dev libhdf5-mpich-dev
#RUN DEBIAN_FRONTEND=noninteractive apt-get install -y clang
#RUN DEBIAN_FRONTEND=noninteractive apt-get install -y texlive
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y ninja-build
RUN git clone https://github.com/GlobalArrays/ga && cd ga && git checkout 2518e23433385bfa3726d507b8cd0d7ed038021b && ./autogen.sh && ./configure --disable-f77 --with-mpi3 && make install && cd .. && rm -rf ga
#RUN DEBIAN_FRONTEND=noninteractive apt-get install -y wget gnupg2 software-properties-common
#RUN wget https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB -O - | apt-key add –
#RUN wget https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB  && apt-key add GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB
#RUN echo "deb https://apt.repos.intel.com/oneapi all main" | tee /etc/apt/sources.list.d/oneAPI.list
#RUN add-apt-repository "deb https://apt.repos.intel.com/oneapi all main"
#RUN DEBIAN_FRONTEND=noninteractive apt install intel-basekit intel-hpckit
#RUN DEBIAN_FRONTEND=noninteractive apt install intel-hpckit
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y wget unzip
RUN wget https://github.com/ninja-build/ninja/releases/download/v1.10.2/ninja-linux.zip && unzip -o -d /usr/bin ninja-linux.zip  && rm ninja-linux.zip
RUN git clone https://github.com/pjknowles/libtensor && cd libtensor && git checkout gmb && mkdir cmake-build && cd cmake-build && cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=ON -DLIBTENSOR_STANDALONE=ON -DUSE_MKL=1 -G Ninja .. && ninja && ninja install && cd ../.. && rm -rf libtensor

RUN apt-get install software-properties-common -y 
RUN add-apt-repository ppa:ubuntu-toolchain-r/test -y 
RUN apt-get update
RUN apt-get install gcc-11 g++-11 gfortran-11 clang-10 -y 
RUN update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-11 60 --slave /usr/bin/g++ g++ /usr/bin/g++-11
RUN update-alternatives --config gcc
