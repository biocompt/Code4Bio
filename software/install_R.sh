#!/bin/sh

NPROC=$(nproc)

# R Installation and Essential Libraries Setup on Ubuntu (22.04 LTS)
# This file provides a step-by-step guide to install R (version 4.2.2) on Ubuntu, along with the essential system libraries required for installing and running R packages.
# It ensures that R and its dependencies are properly configured, enabling seamless installation and usage of additional R packages.

# 1. Install mandatory packages
echo "Installing dependencies..."
sudo apt-get update && sudo apt-get install \
build-essential fort77 xorg-dev liblzma-dev libblas-dev gfortran \
gcc-multilib gobjc++ aptitude libssl-dev libclang-dev \
libreadline-dev libbz2-dev libpcre2-dev \
libcurl4 libcurl4-openssl-dev gdebi-core \
default-jre default-jdk openjdk-8-jdk openjdk-8-jre -y

# 2. Download R-4.2.2
echo "Downloading R-4.2.2..."
wget https://cloud.r-project.org/src/base/R-4/R-4.2.2.tar.gz
tar -xvf R-4.2.2.tar.gz
cd R-4.2.2/

# 3. Install R
echo "Compiling and installing R..."
sudo ./configure CFLAGS='-g -O0' CXXFLAGS='-g -O0' --enable-R-shlib --enable-memory-profiling --without-recommended-packages
make  -j "$NPROC"
make check
sudo make install

# 4. Optional: install Rstudio
echo "Installing RStudio..."
wget https://download1.rstudio.org/electron/jammy/amd64/rstudio-2024.12.1-563-amd64.deb
sudo gdebi rstudio-2024.12.1-563-amd64.deb

echo "Installation completed successfully!"