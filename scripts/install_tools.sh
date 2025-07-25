#!/bin/bash

# Install (bcftools, htslib, trf)

set -e

# Initialize flags
INSTALL_BCFTOOLS=false
INSTALL_HTSLIB=false
INSTALL_TRF=false
INSTALL_PARALLEL=false

# Parse arguments and set flags
for tool in "$@"; do
    case "$tool" in
        bcftools)
            INSTALL_BCFTOOLS=true
            ;;
        htslib)
            INSTALL_HTSLIB=true
            ;;
        trf)
            INSTALL_TRF=true
            ;;
        *)
            echo "Unknown tool: $tool"
            echo "Usage: $0 [bcftools] [htslib] [trf]"
            exit 1
            ;;
    esac
done

# Run actual installation based on flags
if $INSTALL_BCFTOOLS; then
    echo "Installing bcftools..."
    wget -q https://github.com/samtools/bcftools/releases/download/1.21/bcftools-1.21.tar.bz2
    tar -xf bcftools-1.21.tar.bz2
    cd bcftools-1.21/
    ./configure
    make -j8
    cd ..
    rm bcftools-1.21.tar.bz2
    echo "done installing bcftools"
fi

if $INSTALL_HTSLIB; then
    echo "Installing htslib..."
    wget -q https://github.com/samtools/htslib/releases/download/1.21/htslib-1.21.tar.bz2
    tar -xf htslib-1.21.tar.bz2
    cd htslib-1.21/
    ./configure
    make -j8
    cd ..
    rm htslib-1.21.tar.bz2
    echo "done installing htslib"
fi

if $INSTALL_TRF; then
    echo "Installing TRF..."
    wget -q https://github.com/Benson-Genomics-Lab/TRF/releases/download/v4.09.1/trf409.linux64
    chmod +x trf409.linux64
    echo "done installing trf"
fi

if $INSTALL_PARALLEL; then
    echo "Installing GNU Parallel..."
    wget -q https://ftp.gnu.org/gnu/parallel/parallel-latest.tar.bz2
    tar -xf parallel-latest.tar.bz2
    cd parallel-*
    ./configure
    make -j8
    cd ..
    rm parallel-latest.tar.bz2
    echo "done installing GNU Parallel"
fi

# uncomment to install tools required for simulating SVs (samtools, VISOR, sniffles, minimap2, pbsim3, ccs)

# wget https://github.com/samtools/samtools/releases/download/1.21/samtools-1.21.tar.bz2
# tar -xvf samtools-1.21.tar.bz2
# cd samtools-1.21/
# ./configure
# make -j 8
# cd ..
# rm samtools-1.21.tar.bz2

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
