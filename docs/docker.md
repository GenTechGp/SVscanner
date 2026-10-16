# SVscanner container

The container image is the supported way to run SVscanner. It is published to the GitHub
Container Registry on every release tag:

```
ghcr.io/gentechgp/svscanner:X.Y.Z     # a specific release
ghcr.io/gentechgp/svscanner:X.Y       # latest patch of a minor series
ghcr.io/gentechgp/svscanner:latest    # latest release
```

Pin an exact version for anything whose results you intend to keep.

## What is in it

Built on [`dfam/tetools`](https://github.com/Dfam-consortium/TETools), the container image
maintained by the Dfam consortium, which supplies a correctly configured RepeatMasker.

| Component | Version | Source |
|---|---|---|
| RepeatMasker | 4.2.6 | `dfam/tetools:2.01` |
| Dfam database | 4.0: base file + curated consensus partition (FamDB 3.0) | `dfam/tetools:2.01` |
| RMBlast (default search engine) | 2.17.1 | `dfam/tetools:2.01` |
| HMMER | 3.4 | `dfam/tetools:2.01` |
| TRF | 4.09.1 | `dfam/tetools:2.01` |
| bcftools, bgzip, tabix | 1.21 | built from source |
| GNU parallel, Python 3.11 + dependencies | Debian 12 | apt / pip |

The resolved Python package versions for a given image are recorded inside it at
`/opt/svscanner-venv/requirements.lock`.

SVscanner requires RepeatMasker 4.2.6 or later and checks this at startup. Earlier
versions split forward-strand SVA insertions with expanded VNTRs into several fragments
while reporting the reverse-strand copy as one element
([Dfam-consortium/RepeatMasker#386](https://github.com/Dfam-consortium/RepeatMasker/issues/386)).

The image is around 5.5 GB unpacked, nearly all of it the `dfam/tetools` base, which also
carries RepeatModeler, MAFFT, genometools and other tools SVscanner does not use. The base
is Dfam-maintained and correctly configured, which was judged a better trade than
maintaining a smaller RepeatMasker build.

## The Dfam database

The built-in Dfam 4.0 curated consensus partition is all RMBlast needs. It holds the
curated families of every species, so `--species human` resolves 1,387 families (and, say,
`--species "danio rerio"` 1,966) from the image alone, with nothing to download or mount.

HMMER, uncurated families or another Dfam release need more of Dfam (see
[Dfam releases](https://www.dfam.org/releases/)). Mount a directory of FamDB files into the
container and pass it with one of:

- `--dfam_extra DIR` to **add** partitions of the same release (Dfam 4.0) to the built-in
  database, e.g. `dfam40.curated.hmm.*.h5` for `--engine hmmer`, or
  `dfam40.uncurated.consensus.*.h5` together with `--uncurated`;
- `--dfam_db DIR` to **replace** the built-in database with a complete one (its base file
  `dfamNN.0.h5` plus partitions), e.g. a later Dfam release.

Details and the checks SVscanner makes are in [Dfam databases](Commands.md#dfam-databases).

Each run spends about 20 seconds preparing RepeatMasker's search library before
annotating, reported in the log as `Library preparation took N seconds`. The library is
built inside the run's output directory (`REPEATMASKER_CACHE_DIR`) and removed at the end,
so nothing is written into the image or your home directory and concurrent runs never
share it.

## Docker

```
docker run --rm \
  -v /path/to/reference:/ref:ro \
  -v "$PWD":/data \
  ghcr.io/gentechgp/svscanner:X.Y.Z \
  svscanner --vcf /data/input.vcf.gz --ref /ref/hg38.fa --out /data/svscanner_out \
            --nthread 16
```

The image's working directory is `/data`. Outputs are written as the container's user:
root under standard Docker, which is your own user under rootless Docker. With standard
Docker, add `--user "$(id -u):$(id -g)"` to have the outputs owned by you; under rootless
Docker do not, as it makes the mounted directories unwritable.

## Singularity / Apptainer

```
singularity pull svscanner_X.Y.Z.sif docker://ghcr.io/gentechgp/svscanner:X.Y.Z

singularity exec \
  --bind /path/to/reference:/ref:ro \
  svscanner_X.Y.Z.sif \
  svscanner --vcf input.vcf.gz --ref /ref/hg38.fa --out svscanner_out --nthread 16
```

The container filesystem is read only under Singularity. Nothing in SVscanner writes back
into the image, so this needs no special handling. On a cluster, pull the `.sif` on a node
with internet access and keep it on shared storage; see the README for a PBS job script.

## Nextflow

SVscanner does not ship a Nextflow pipeline, but the image works as a process container
as is: `svscanner` is on the image's `PATH` and nothing depends on an `ENTRYPOINT`, which
Nextflow overrides. A minimal pipeline:

`main.nf`

```groovy
nextflow.enable.dsl = 2

params.vcf        = null          // SV VCF
params.ref        = null          // reference FASTA; its .fai must sit next to it
params.species    = 'human'
params.dfam_extra = null          // optional: directory of extra Dfam partitions
params.outdir     = 'results'

process SVSCANNER {
    tag { vcf.simpleName }
    publishDir params.outdir, mode: 'copy'

    input:
      path vcf
      tuple path(ref), path(ref_fai)
      path dfam_extra             // an empty list when not given

    output:
      path "${vcf.simpleName}_svscanner"

    script:
      def dfamArg = dfam_extra ? "--dfam_extra ${dfam_extra}" : ""
      """
      svscanner --vcf ${vcf} --ref ${ref} \\
                --out ${vcf.simpleName}_svscanner \\
                --species '${params.species}' \\
                --nthread ${task.cpus} ${dfamArg}
      """
}

workflow {
    if (!params.vcf || !params.ref) error "Set --vcf and --ref"
    SVSCANNER(
        file(params.vcf),
        tuple(file(params.ref), file("${params.ref}.fai")),
        params.dfam_extra ? file(params.dfam_extra) : []
    )
}
```

`nextflow.config`

```groovy
process {
    withName: 'SVSCANNER' {
        container = 'ghcr.io/gentechgp/svscanner:X.Y.Z'   // or a local .sif path
        cpus   = 16
        memory = '64GB'
        time   = '12h'
    }
}

profiles {
    docker      { docker.enabled = true }
    singularity { singularity.enabled = true; singularity.autoMounts = true }
}
```

Run it:

```
nextflow run main.nf -profile singularity --vcf input.vcf.gz --ref /path/hg38.fa
```

On a cluster whose compute nodes have no internet access, point `container` at a `.sif`
you pulled beforehand instead of the registry reference.

Two things are worth getting right:

- **Pass `--nthread ${task.cpus}`.** SVscanner otherwise sizes itself from `$PBS_NCPUS`,
  or from `nproc` when that is unset. `nproc` honours a cpuset but not a CPU *quota*, so
  inside a container limited with `--cpus` it can report the whole host and oversubscribe
  badly. There is no `$PBS_NCPUS` inside the container even on a PBS cluster, because the
  process is not a direct child of the job script.
- **Stage everything SVscanner reads as an input `path`**, including the reference index
  and any `--dfam_extra` / `--dfam_db` directory. Nextflow then mounts them under Docker
  and Singularity alike; paths that are not staged are invisible inside the container.

## Building locally

```
git archive HEAD | docker build -t svscanner:dev-$(git rev-parse --short HEAD) -
```

`git archive` builds from the committed tree only, so untracked files cannot leak into the
image. To get a `.sif` from it (with rootless Docker, `apptainer` needs `DOCKER_HOST`):

```
export DOCKER_HOST=unix:///run/user/$(id -u)/docker.sock
apptainer build svscanner_dev.sif docker-daemon://svscanner:dev-$(git rev-parse --short HEAD)
```

The build runs its own checks: Python dependencies, every required binary, RepeatMasker
4.2.6 or later, and that `famdb.py` can start and read the built-in database. A build that
succeeds is one where the tooling is wired up correctly. CI additionally runs RepeatMasker
on an SVA in both orientations (see
[`docker-publish.yml`](../.github/workflows/docker-publish.yml)). Neither proves the
annotation of a real dataset is correct.

The base image tag is pinned in the [Dockerfile](../Dockerfile). Read the comment there
before changing it.
