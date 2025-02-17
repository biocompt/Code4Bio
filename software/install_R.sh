#!/bin/sh

# R Installation and Essential Libraries Setup on Ubuntu
# This file provides a step-by-step guide to install R (version 4.2.2) on Ubuntu, along with the essential system libraries required for installing and running R packages.
# It ensures that R and its dependencies are properly configured, enabling seamless installation and usage of additional R packages.

# 1. Install mandatory packages
sudo apt-get update && sudo apt-get install \
build-essential fort77 xorg-dev liblzma-dev libblas-dev gfortran \
gcc-multilib gobjc++ aptitude \
libreadline-dev libbz2-dev libpcre2-dev \
libcurl4 libcurl4-openssl-dev \
default-jre default-jdk openjdk-8-jdk openjdk-8-jre -y

# 2. Download R-4.2.2
wget https://cloud.r-project.org/src/base/R-4/R-4.2.2.tar.gz
tar -xvf R-4.2.2.tar.gz
cd R-4.2.2/

# 3. Install R
sudo ./configure CFLAGS='-g -O0' CXXFLAGS='-g -O0' --enable-R-shlib --enable-memory-profiling --without-recommended-packages
make  -j 16
make check
sudo make install

# 4. Optional: install Rstudio
# Ubuntu 22: https://download1.rstudio.org/electron/jammy/amd64/rstudio-2024.12.1-563-amd64.deb. Then install with Software
