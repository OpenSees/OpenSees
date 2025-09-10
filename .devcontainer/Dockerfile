FROM quay.io/pypa/manylinux_2_28_x86_64

RUN yum install -y \
    openmpi-devel \
    && yum clean all \
    && rm -rf /var/cache/yum

# Add OpenMPI to path
ENV PATH="/usr/lib64/openmpi/bin:$PATH"
ENV LD_LIBRARY_PATH="/usr/lib64/openmpi/lib:$LD_LIBRARY_PATH"

# Get the latest OpenBLAS and build it (v0.3.31.dev as of 2026-02-25 to get BUILD_DOUBLE make fixes)
RUN git clone --depth 1 --revision 18638c70eff1f0d08e2833b2724deaa128d6a334 https://github.com/OpenMathLib/OpenBLAS.git \
    && cd OpenBLAS \
    && make_flags=( \
        DYNAMIC_ARCH=1 \
        TARGET=NEHALEM \
        DYNAMIC_LIST="NEHALEM SANDYBRIDGE HASWELL SKYLAKEX" \
        BUILD_DOUBLE=1 \
        NO_AFFINITY=1 \
        USE_THREAD=1 \
        USE_OPENMP=1 \
        FC=gfortran \
        GEMM_MULTITHREAD_THRESHOLD=50 \
        NUM_THREADS=128 \
        PREFIX=/usr \
        OPENBLAS_INCLUDE_DIR=/usr/include \
        OPENBLAS_LIB_DIR=/usr/lib64 \
        -j$(nproc) \
    ) \
    && make "${make_flags[@]}" \
    && make "${make_flags[@]}" install \
    && cd .. \
    && rm -rf OpenBLAS

ARG USERNAME=u
ARG USER_UID=1000
ARG USER_GID=$USER_UID

RUN groupadd --gid $USER_GID $USERNAME \
    && useradd --uid $USER_UID --gid $USER_GID -m $USERNAME \
    && yum install -y sudo \
    && yum clean all \
    && rm -rf /var/cache/yum \
    && echo $USERNAME ALL=\(root\) NOPASSWD:ALL > /etc/sudoers.d/$USERNAME \
    && chmod 0440 /etc/sudoers.d/$USERNAME

USER $USERNAME