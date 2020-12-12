FROM ubuntu:16.04

COPY . /usr/sian/

RUN  apt-get update -y
RUN  apt-get install curl wget -y 
RUN  apt-get install -y bzip2 gcc cmake gfortran libpcre3-dev build-essential

RUN cd /usr/local/src && \
    wget https://julialang-s3.julialang.org/bin/linux/x64/1.5/julia-1.5.3-linux-x86_64.tar.gz && \
    tar zxvf julia-1.5.3-linux-x86_64.tar.gz

RUN ln -s /usr/local/src/julia-1.5.3/bin/julia /usr/local/bin/julia

ENV JULIA_PKGDIR /root/.julia/
RUN echo "using Pkg\n\
    metadata_packages = [ \"Oscar\", \"LinearAlgebra\", \"Singular\", \"GroebnerBasis\", \"MacroTools\", \"OrderedCollections\"]\n\
    for\n\
    package = metadata_packages \n\
    Pkg.add(package)\n\
    end;\n\ 
    Pkg.resolve()\n" >> /tmp/package_installs.jl
RUN cat /tmp/package_installs.jl
RUN julia /tmp/package_installs.jl
WORKDIR /usr/sian
CMD ["/bin/bash"]