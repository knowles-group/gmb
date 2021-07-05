FROM intel/oneapi-hpckit

RUN apt-get update 
RUN apt-get upgrade -y
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y cmake git doxygen graphviz rsync curl libeigen3-dev wget unzip software-properties-common
RUN wget https://github.com/ninja-build/ninja/releases/download/v1.10.2/ninja-linux.zip && unzip -o -d /usr/bin ninja-linux.zip  && rm ninja-linux.zip

RUN add-apt-repository ppa:ubuntu-toolchain-r/test -y
RUN apt-get update
RUN apt-get install gcc-11 g++-11 gfortran-11 clang-10 -y
RUN update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-11 60 --slave /usr/bin/g++ g++ /usr/bin/g++-11
RUN update-alternatives --config gcc

RUN wget https://github.com/GlobalArrays/ga/releases/download/v5.8/ga-5.8.tar.gz && tar xzf ga-5.8.tar.gz && cd ga-5.8 && ./configure --disable-f77 --with-mpi3 && make install && cd .. && rm -rf ga-5.8*
RUN git clone https://github.com/pjknowles/libtensor && cd libtensor && git checkout gmb && mkdir cmake-build && cd cmake-build && cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=ON -DLIBTENSOR_STANDALONE=ON -DUSE_MKL=1 -G Ninja .. && ninja && ninja install && cd ../.. && rm -rf libtensor
RUN rm -rf /opt/intel/oneapi/{ipp,intelpython,dnnl,conda_channel,itac,vtune,dal,advisor,clck,inspector,ccl}
