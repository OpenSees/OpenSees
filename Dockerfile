# Dockerfile to build OpenSees
#   .. utilizes ubuntu:20.04 LTS as base
#   .. it will build sequential version and place in /usr/local/bin

# written: fmk

FROM ubuntu:20.04

SHELL ["/bin/bash", "-c"]

WORKDIR /opensees

ARG versionOpenSees=v3.6.0

RUN cp /etc/apt/sources.list /etc/apt/sources.list~ \
    && sed -Ei 's/^# deb-src /deb-src /' /etc/apt/sources.list \
    && apt-get update \
    && DEBIAN_FRONTEND=noninteractive TZ=Etc/UTC apt-get -y install tzdata git \
    && apt-get install -y sudo \
    && sudo apt install -y cmake gcc g++ gfortran liblapack-dev git python3-pip \
    && pip3 install conan==1.60.1 \
    && git clone --depth 1 --branch hdf5-1_12_2 https://github.com/HDFGroup/hdf5.git \
    && cd hdf5 \
    && ./configure --prefix=/usr/local/hdf5 \
    && make \
    && cd .. \
    && git clone -b $versionOpenSees --single-branch https://github.com/OpenSees/OpenSees.git \
    && cd OpenSees \
    && mkdir build \
    && cd build \
    && conan install .. --build missing \
    && cmake .. \
    && cmake --build . --config Release \
    && cmake --install . \
    && cp -r ./lib/tcl8.6 /usr/local/lib \
    && cd ../.. \
    && rm -fr OpenSees \
    && rm -fr hdf5

  

    