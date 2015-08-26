# HiFive latest version docker container build file

# Load base container ubuntu version 14.04

FROM ubuntu:14.04

# Load python2.7, samtools, and gfortran

RUN apt-get update
RUN apt-get install -y python2.7-dev samtools gfortran

# Load libraries hdf5, atlas, 

RUN apt-get install -y libhdf5-dev libatlas-base-dev

# Load mpi and supporting library

# RUN apt-get install openmpi-bin openmpi-common openssh-client openssh-server libopenmpi1.6 libopenmpi-dev -y

# Get pip

RUN apt-get install -y python-pip
RUN pip install -U pip
RUN pip install -U setuptools

# Get python packages: numpy, scipy, pysam, cython, mpi4py, and h5py

# RUN pip install numpy scipy pysam cython mpi4py
RUN pip install numpy scipy pysam cython
RUN pip install h5py

# Get hifive

RUN apt-get install -y unzip wget
RUN mkdir -p tmp && cd tmp && wget https://github.com/bxlab/hifive/archive/master.zip && unzip master.zip && cd hifive-master && python setup.py install

MAINTAINER Michael Sauria <mike.sauria@jhu.edu>
