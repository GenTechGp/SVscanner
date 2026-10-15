# SVscanner container
#
# Built on the Dfam TE Tools image, which supplies a configured RepeatMasker together
# with TRF, RMBlast, HMMER and FamDB. The tag is pinned deliberately:
#
#   dfam/tetools:2.01 -> RepeatMasker 4.2.6, TRF 4.09.1, RMBlast 2.17.1, HMMER 3.4,
#                        FamDB 3.0.0 with Dfam 4.0 (base + curated consensus partitions)
#
# SVscanner needs RepeatMasker 4.2.6+ (run_workflow.sh checks): it fixes forward-strand
# SVAs being split into fragments, and adds the FAMDB_DATA_DIR / REPEATMASKER_CACHE_DIR
# settings the workflow relies on.
#
# The curated consensus partition is all RMBlast needs: `-species human` resolves 1,387
# curated families from the image alone. Other engines or species add partitions at
# runtime with --dfam_extra, or replace the database with --dfam_db.
FROM dfam/tetools:2.01

LABEL org.opencontainers.image.title="SVscanner" \
      org.opencontainers.image.description="Annotate tandem repeats and mobile elements within structural variants using TRF and RepeatMasker" \
      org.opencontainers.image.source="https://github.com/GenTechGp/SVscanner" \
      org.opencontainers.image.base.name="docker.io/dfam/tetools:2.01" \
      org.opencontainers.image.licenses="MIT"

ARG HTSLIB_VERSION=1.21
ARG BCFTOOLS_VERSION=1.21

# GNU parallel drives the TRF and RepeatMasker fan-out; python3-venv is needed
# because Debian 12 marks its system Python as externally managed (PEP 668).
RUN apt-get -y update && apt-get -y install --no-install-recommends \
        parallel \
        python3-venv \
        libbz2-dev \
        liblzma-dev \
        libcurl4-openssl-dev \
        libssl-dev \
    && rm -rf /var/lib/apt/lists/*

# bcftools, bgzip and tabix. Built from source rather than taken from Debian 12
# (which has 1.16). htslib's lib directory has to go on the loader path: bcftools
# links libhts dynamically without an rpath, so without ldconfig it builds fine and
# then fails to start.
RUN cd /tmp \
    && curl -fsSL -o htslib.tar.bz2 "https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2" \
    && tar -xf htslib.tar.bz2 \
    && cd "htslib-${HTSLIB_VERSION}" \
    && ./configure --prefix=/opt/htslib \
    && make -j"$(nproc)" && make install \
    && cd /tmp && rm -rf htslib.tar.bz2 "htslib-${HTSLIB_VERSION}" \
    && echo /opt/htslib/lib > /etc/ld.so.conf.d/htslib.conf && ldconfig

RUN cd /tmp \
    && curl -fsSL -o bcftools.tar.bz2 "https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2" \
    && tar -xf bcftools.tar.bz2 \
    && cd "bcftools-${BCFTOOLS_VERSION}" \
    && ./configure --prefix=/opt/bcftools --with-htslib=/opt/htslib \
    && make -j"$(nproc)" && make install \
    && cd /tmp && rm -rf bcftools.tar.bz2 "bcftools-${BCFTOOLS_VERSION}"

# Python dependencies, in their own layer so source changes do not reinstall them.
# The resolved versions are recorded in the image: requirements.txt is deliberately
# loose, so this is the only record of what a given image tag actually contains.
COPY requirements.txt /opt/svscanner/requirements.txt
RUN python3 -m venv /opt/svscanner-venv \
    && /opt/svscanner-venv/bin/pip install --no-cache-dir --upgrade pip \
    && /opt/svscanner-venv/bin/pip install --no-cache-dir -r /opt/svscanner/requirements.txt \
    && /opt/svscanner-venv/bin/pip freeze > /opt/svscanner-venv/requirements.lock \
    && cat /opt/svscanner-venv/requirements.lock

COPY . /opt/svscanner

# `svscanner` on $PATH is the documented entry point. Nextflow overrides a container's
# ENTRYPOINT, so nothing may depend on one: every path below is set through ENV instead.
RUN printf '#!/bin/bash\nexec /opt/svscanner/scripts/run_workflow.sh "$@"\n' > /usr/local/bin/svscanner \
    && chmod +x /usr/local/bin/svscanner

ENV PATH=/opt/svscanner-venv/bin:/opt/bcftools/bin:/opt/htslib/bin:$PATH \
    SVSCANNER_TESTDATA=/opt/svscanner/test \
    MPLCONFIGDIR=/tmp/matplotlib

# Fail the build rather than the pipeline. Covers the things most likely to break on a
# base image bump: a missing tool, a RepeatMasker older than run_workflow.sh accepts, and
# famdb.py unable to start. famdb.py is executed directly, the way RepeatMasker runs it,
# so its own interpreter line decides which Python (and h5py) is used.
#
# Two things this has to get right, both learned the hard way:
#   - one `command -v` per tool. It accepts several names but reports success if any one
#     of them resolves, which waves through a half-broken image.
#   - no pipes. This runs under dash, which has no `set -o pipefail`, so `cmd | head -1`
#     reports head's exit status and hides a tool that cannot even start.
RUN set -eu; \
    for tool in trf RepeatMasker parallel bcftools bgzip tabix python3; do \
        command -v "$tool" > /dev/null || { echo "missing from PATH: $tool" >&2; exit 1; }; \
    done; \
    /opt/svscanner/scripts/check_required_python.sh; \
    bcftools --version; \
    bgzip --version; \
    tabix --version; \
    RepeatMasker -v; \
    rm_version=$(RepeatMasker -v | awk '/RepeatMasker version/ { print $3 }'); \
    [ "$(printf '4.2.6\n%s\n' "${rm_version}" | sort -V | head -n1)" = "4.2.6" ] \
        || { echo "RepeatMasker '${rm_version}' is older than 4.2.6" >&2; exit 1; }; \
    /opt/FamDB/famdb.py info; \
    svscanner --version

WORKDIR /data
