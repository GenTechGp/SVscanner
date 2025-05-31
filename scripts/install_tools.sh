#!/bin/bash

wget https://github.com/samtools/samtools/releases/download/1.21/samtools-1.21.tar.bz2
tar -xvf samtools-1.21.tar.bz2
cd samtools-1.21/
./configure
make -j 8
cd ..
rm samtools-1.21.tar.bz2

wget https://github.com/samtools/bcftools/releases/download/1.21/bcftools-1.21.tar.bz2
tar -xvf bcftools-1.21.tar.bz2
cd bcftools-1.21/
./configure
make -j 8
cd ..
rm bcftools-1.21.tar.bz2

wget https://github.com/samtools/htslib/releases/download/1.21/htslib-1.21.tar.bz2
tar -xvf htslib-1.21.tar.bz2
cd htslib-1.21/
./configure
make -j 8
cd ..
rm htslib-1.21.tar.bz2

wget https://github.com/Benson-Genomics-Lab/TRF/releases/download/v4.09.1/trf409.linux64
chmod +x trf409.linux64

# uncomment to install tools required for simulating SVs
# wget https://github.com/davidebolo1993/VISOR/archive/refs/tags/v1.1.2.1.tar.gz
# tar -xvf v1.1.2.1.tar.gz
# rm v1.1.2.1.tar.gz
# cd VISOR-1.1.2.1/
# python3.10 -m venv visor
# source visor/bin/activate
# pip install -r requirements.txt
# python setup.py install
# VISOR --help
# deactivate
# cd ..

# python3.10 -m venv sniffles_260
# source sniffles_260/bin/activate
# pip install sniffles
# sniffles --version
# deactivate

# wget https://github.com/lh3/minimap2/releases/download/v2.29/minimap2-2.29_x64-linux.tar.bz2
# tar -xvf minimap2-2.29_x64-linux.tar.bz2
# rm minimap2-2.29_x64-linux.tar.bz2

# wget https://github.com/yukiteruono/pbsim3/archive/refs/tags/v3.0.5.tar.gz
# tar -xvf v3.0.5.tar.gz
# rm v3.0.5.tar.gz
# cd pbsim3-3.0.5/
# ./configure
# make -j 8
# cd ..

# mkdir ccs_v6.4.0 && cd ccs_v6.4.0
# wget https://github.com/PacificBiosciences/ccs/releases/download/v6.4.0/ccs.tar.gz
# tar -xvf ccs.tar.gz
# rm ccs.tar.gz
# cd ..
