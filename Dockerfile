FROM ubuntu:16.04

COPY examples/ /usr/sian/examples/
COPY without-macros/ /usr/sian/without-macros/
COPY IdentifiabilityODE.jl /usr/sian
COPY SIAN /usr/sian/SIAN
ADD package_installs.jl /tmp/package_installs.jl

RUN  apt-get update -y
RUN  apt-get install curl wget -y 
# git software-properties-common libcairo2 libpango1.0-0 -y
# RUN  apt-get install -y libpcre3-dev build-essential
# RUN  apt-get install -y gettext hdf5-tools
# RUN  apt-get install -y gfortran python
# RUN  apt-get install -y m4 cmake libssl-dev 
RUN  apt-get install -y bzip2 gcc  cmake gfortran libpcre3-dev build-essential

RUN cd /usr/local/src && \
    wget https://julialang-s3.julialang.org/bin/linux/x64/1.5/julia-1.5.3-linux-x86_64.tar.gz && \
    tar zxvf julia-1.5.3-linux-x86_64.tar.gz

RUN ln -s /usr/local/src/julia-1.5.3/bin/julia /usr/local/bin/julia

ENV JULIA_PKGDIR /root/.julia/

RUN julia /tmp/package_installs.jl

# RUN julia -e "using Oscar;" && \
#     julia -e "using LinearAlgebra;" && \
#     julia -e "using Singular;" && \
#     julia -e "using GroebnerBasis;" && \
#     julia -e "using MacroTools;" && \
#     julia -e "using OrderedCollections;"

WORKDIR /usr/sian
CMD ["/bin/bash"]