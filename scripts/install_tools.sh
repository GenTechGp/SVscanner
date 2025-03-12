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