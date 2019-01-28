################## BASE IMAGE ###################### 
FROM continuumio/miniconda3:latest  

################## METADATA ###################### 
LABEL about.summary="WGS-EQA: nextflow based pipeline for external QA for bacterial WGS"

################## MAINTAINER ###################### 
MAINTAINER David Eyre <david.eyre@bdi.ox.ac.uk>  

################## INSTALLATION ######################  

RUN apt-get clean all && \
 apt-get update && \
 apt-get upgrade -y && \
 apt-get install -y procps && \
 apt-get clean && \
 apt-get purge && \
 rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge

RUN conda install wgsim=1.0 \
				  numpy=1.15.4 \
				  biopython=1.72

WORKDIR /data/

