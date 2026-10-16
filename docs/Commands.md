# Commands
1. [run_workflow.sh](#run_workflowsh)
2. [extract_sv.py](#extract_svpy)
3. [repeat_annotatoin.py](#repeat_annotatoinpy)
4. [generate_plot.py](#generate_plotpy)
5. [simulate_sv.py](#simulate_svpy)
6. [check_vcf_ids.py](#check_vcf_idspy)
7. [check_ref.py](#check_refpy)

## run_workflow.sh

---

**Required Arguments**

| Argument       | Type   | Description                                                        |
|----------------|--------|--------------------------------------------------------------------|
| `--out`        | `str`  | **Path to the output directory**.                                 |
| `--vcf`        | `str`  | **Path to the structural variant (SV) VCF file**.                 |
| `--ref`        | `str`  | **Path to the reference FASTA file**.                             |

---

**Optional Arguments**

| Argument                        | Type           | Default                 | Description                                                                 |
|---------------------------------|----------------|-------------------------|-----------------------------------------------------------------------------|
| `--prefix`                      | `str`          | `None`                  | Prefix for output files.                                                    |
| `--str_bed`                     | `str`          | `None`              | Path to BED file containing STR (short tandem repeat) elements.            |
| `--species`                     | `str`          | `human`              | Species name used by RepeatMasker.                                         |
| `--engine`                      | `str`          | RepeatMasker's default (RMBlast in the container) | RepeatMasker search engine, e.g. `rmblast`, `hmmer`. |
| `--dfam_db`                     | `str`          | `None`                  | A complete Dfam database (base file `dfamNN.0.h5` plus partitions) to use instead of the built-in one. See [Dfam databases](#dfam-databases). |
| `--dfam_extra`                  | `str`          | `None`                  | Dfam partitions to add to the built-in database (same Dfam release, no base file). See [Dfam databases](#dfam-databases). |
| `--uncurated`                   | flag           |                         | Use uncurated Dfam families as well as curated ones (RepeatMasker `-uncurated`). Needs uncurated partitions covering `--species`. See [Dfam databases](#dfam-databases). |
| `--rm_gc`                       | `int`          | not set (RepeatMasker chooses) | Background GC% for RepeatMasker's scoring matrices (RepeatMasker `-gc`, 1-100). `--rm_gc 43` makes annotations independent of how SVs are grouped into files; see the README Notes. |
| `--min_sv_coverage`             | `float`        | `0.05`      | Minimum intersection between a repeat element and SV.                      |
| `--min_class_sv_coverage`       | `float`        | `0.25`| Minimum class-level SV coverage to be considered repetitive.               |
| `--min_total_sv_coverage`       | `float`        | `0.75`| Minimum total SV coverage by repeats to be considered repetitive.          |
| `--max_trf_overlap`             | `float`        | `0.1`      | Maximum TRF element overlap to be considered non-overlapping.              |
| `--interval`                    | `int`          | `0.05`             | Interval value for binning or windowing.                                   |
| `--diagram_len`                 | `int`          | `100`          | Length of the diagram generated for visualization.                         |
| `--no_rm_chain`                 | flag           |                | Do not chain RepeatMasker alignments per element: RepeatMasker runs without `-a`, no `rm_chain.tab` is written and the diagram has no `RM chain` block. Default: chain. |
| `--rm_chain_jump_cost`          | `float`        | `30`           | Cost, in matched bases, of each jump between RepeatMasker alignment blocks when they are chained into one alignment per element (the `RM chain` diagram block, see [rm_align_chain.py](#rm_align_chainpy)). `--resume` reuses the existing `rm_chain.tab`, so a different value needs a full run. |
| `--nsplit_files`                | `int`          | `500`         | Number of files to split sequences into.                                   |
| `--keep_tmp_files`             | flag           |                         | Keep intermediate files (default: delete after run).                       |
| `--overwrite`                  | flag           |                         | Overwrite existing output files if present.                                |
| `--populate_ids`               | flag           |                         | Fill missing and duplicated IDs in the VCF ID column, in a copy of the VCF. Without it the run exits when the ID column is unusable. See [VCF ID column](#vcf-id-column). |
| `--skip_ref_check`             | flag           |                         | Do not check that `--ref` is the reference the VCF was called against. See [Reference check](#reference-check). |
| `--nthread`                    | `int`          | all available threads   | Number of threads to use.                                                  |
| `--njob`                       | `int`          | `48`             | Number of parallel jobs for RepeatMasker.                                  |
| `--help`                       | flag           |                         | Show help message and exit.                                                |
| `--version`                    | flag           |                         | Show version information and exit.                                         |

---

### VCF ID column

SVscanner adds its annotations back onto the VCF by the ID column, so every record needs an
ID that is present (not `.`) and unique. Before anything else runs, `run_workflow.sh` checks
this with [check_vcf_ids.py](#check_vcf_idspy):

- **IDs fine:** the run continues unchanged.
- **Missing or duplicated IDs, no `--populate_ids`:** the run exits with an error and creates
  no output directory.
- **Missing or duplicated IDs, with `--populate_ids`:** the input VCF is never edited. A
  repaired copy is written to `<out>/<prefix>ids_populated.vcf.gz`, together with
  `<out>/<prefix>ids_populated.map.tsv` (old ID to new ID), and that copy is annotated.
  Existing unique IDs are kept; the first occurrence of a duplicated ID is kept. Every
  missing ID and every later duplicate becomes `SVSCANNER_<CHROM>_<POS>_<SVTYPE>_<n>`, where
  `<n>` is the 1-based record number. Renamed duplicates are reported as a warning. The
  annotated VCF therefore carries these IDs.
- **Not repairable, with or without the flag:** an ID containing whitespace, or missing or
  duplicated IDs in a VCF whose records carry `MATEID`. Rewriting IDs would leave the
  `MATEID` references pointing at the wrong records, so fix the IDs at the source.

`--resume` only checks. If the previous run used `--populate_ids`, the repaired copy in the
output directory is used automatically (it is what the saved `.tab` files were built from,
so it is kept after the run); otherwise the VCF given in `--vcf` is checked and the run exits
if its IDs are unusable.

After annotation, `annotate_vcf.py` reports how many VCF records matched an annotation and
warns when none did.

### Reference check

`--ref` must be the reference the VCF was called against. A different build does not make
the run fail: the DEL, DUP, INV and BND sequences, and the flanks of every SV, are fetched
from it and annotated as if they were right (typically many DELs come out `NON_REPETITIVE`).
Right after the ID check, and before the output directory is created, `run_workflow.sh`
checks the reference with [check_ref.py](#check_refpy). The run exits with an error when:

- a contig the records use (CHROM, or a BND mate in ALT) is not in the reference;
- a `##contig` length in the VCF header differs from the reference length for that contig;
- more than 5% of a random sample of REF alleles differ from the reference at their position.

`--skip_ref_check` turns the check off. `--resume` does not run it (no sequence is extracted).

### Dfam databases

RepeatMasker reads its Dfam families through FamDB. A Dfam release (4.0 onwards, FamDB
3.0 format) is a set of `.h5` files:

- one **base file**, `dfamNN.0.h5` (e.g. `dfam40.0.h5`): family names and taxonomy, always
  required;
- **partitions** by component and taxon, e.g. `dfam40.curated.consensus.0.h5` (the curated
  consensus sequences RMBlast searches), `dfam40.curated.hmm.*.h5` (profile HMMs for
  HMMER), `dfam40.uncurated.*.h5` (uncurated families, mainly for species with little
  curated data).

Partitions only work with the base file of their own release. Downloads are on the
[Dfam releases](https://www.dfam.org/releases/) page.

The container has the Dfam 4.0 base file and curated consensus partition built in, which
covers `--species human` with RMBlast. Two options change that, and they cannot be
combined:

| | `--dfam_extra DIR` | `--dfam_db DIR` |
|---|---|---|
| Effect | Adds `DIR`'s partitions to the built-in database | Uses `DIR` instead of the built-in database |
| `DIR` must contain | Partitions of the built-in release only (`dfam40.*`), no base file | Exactly one base file, plus the partitions you need |
| Typical use | `--engine hmmer` (add `dfam40.curated.hmm.*`), uncurated clades | A different Dfam release |

```
# Another species with RMBlast: nothing to add, the built-in database has it
svscanner --species "danio rerio" --vcf ... --ref ... --out ...

# HMMER for human: add the curated HMM partitions to the built-in Dfam 4.0
svscanner --engine hmmer --dfam_extra /data/dfam40_hmm --vcf ... --ref ... --out ...

# A complete database of another release
svscanner --dfam_db /data/dfam41 --species "mus musculus" --vcf ... --ref ... --out ...
```

The built-in curated consensus partition is not limited to human: it is Dfam 4.0's single
root-level partition, holding all 30,659 curated families of every species. With RMBlast,
`--species` alone is therefore enough for any species Dfam has curated families for; the
run log shows how many were found (for `danio rerio`: `249/249 families in ancestor taxa;
1717/1717 lineage-specific families`). HMMER is
different, because its curated HMM partitions are split by taxon and none are built in.

#### Uncurated families

By default RepeatMasker uses curated Dfam families only. `--uncurated` adds Dfam's
uncurated families too, which matters mainly for species with few curated ones. Not every
species has any: in Dfam 4.0 human and `danio rerio` have none, `mus musculus` has 4,692
and `drosophila melanogaster` 599. To check a species before downloading anything (the
base file knows the totals even without the partitions):

```
famdb.py lineage "mus musculus" --uncurated --ancestors --descendants --model consensus -f totals
```

Uncurated families are in their own partitions, none of which are built in, so add the
ones covering your species:

| Search engine | Uncurated partitions to add | Size (Dfam 4.0, compressed) |
|---|---|---|
| RMBlast | `dfam40.uncurated.consensus.0.h5` and `.1.h5` | 0.3 MB and 4.1 GB |
| HMMER | `dfam40.uncurated.hmm.N.h5` for your species' clades, plus the curated HMMs | 1–2.5 GB each |

```
# Mouse, RMBlast, curated + uncurated
svscanner --species "mus musculus" --uncurated --dfam_extra /data/dfam40_uncurated \
          --vcf ... --ref ... --out ...
```

The run log says which families were used (`Including curated and uncurated families:`).
In two cases RepeatMasker would carry on with less than `--uncurated` asked for, and
SVscanner stops the run instead:

- Dfam has no uncurated families for the species (RepeatMasker notes this and runs
  curated-only). Drop `--uncurated`; no partition can help.
- The species has uncurated families but their partitions are not installed (RepeatMasker
  warns that partitions are missing and runs with what is present). Add them with
  `--dfam_extra`. The same check applies without `--uncurated`, e.g. `--engine hmmer`
  with only some of the curated HMM partitions.

With `--dfam_extra`, SVscanner links the built-in files and yours into one directory under
the output directory (symlinks only, removed at the end of the run); a file of yours with
the same name as a built-in one takes its place. With `--dfam_db`, RepeatMasker reads
`DIR` directly. Either way RepeatMasker is pointed at the result through `FAMDB_DATA_DIR`.

Before anything else runs, SVscanner stops with an explanation if:

- `--dfam_extra` contains a base file (use `--dfam_db`), or a partition of a different
  release than the built-in one;
- `--dfam_db` has no base file (use `--dfam_extra`), or more than one;
- FamDB does not actually load one of your files. FamDB skips files it cannot use (for
  example Dfam 3.9 and earlier, which use the older FamDB 2.0 format) without an error, so
  RepeatMasker would otherwise annotate against less than you asked for.

The run log lists every FamDB file in use. As with the built-in database, the species
search library is built per run under the output directory.

## extract_sv.py

---

**Required Arguments**

| Argument              | Type   | Description                                                                 |
|-----------------------|--------|-----------------------------------------------------------------------------|
| `-v`, `--vcf`         | `str`  | **Path to the input VCF file** (can be compressed or uncompressed).        |
| `-r`, `--ref`         | `str`  | **Path to the reference FASTA file**.                                      |
| `-o`, `--out`         | `str`  | **Output directory** to write the resulting FASTA and other files.         |
| `-i`, `--info`        | `str`  | **Path to write the info file** describing processed SVs.                  |

---

**Optional Arguments**

| Argument               | Type           | Default | Description                                                                                       |
|------------------------|----------------|---------|---------------------------------------------------------------------------------------------------|
| `--min`                | `positive int` | `50`    | Minimum SV length to include.                                                                     |
| `--max`                | `positive int` | `50000` | Maximum SV length to include.                                                                     |
| `--flen`               | `positive int` | `2000`  | Max detectable period size supported by TRF (used to determine flanking sequence length).         |
| `--ffac`               | `positive int` | `10`    | Multiplication factor for `SVLEN` to determine flanking sequence length.                          |
| `-n`                   | `positive int` | `1`     | Number of output FASTA files to split sequences evenly into.                                      |
| `--debug`              | flag           |         | Enable debug mode. Prints extra information for troubleshooting.                                  |
| `--warning_count`      | `positive int` | `10`    | Maximum number of warnings to show for each type of warning.                                      |
| `-h`, `--help`         | flag           |         | Show help message and exit.                                                                       |

---

*The `flen` and `ffac` arguments control how much **flanking sequence** is extracted around each SV.

## repeat_annotatoin.py

**Required Arguments**

| Argument      | Type | Required | Description |
|---------------|------|----------|-------------|
| `-v`, `--vcf` | str  | Yes      | Path to the input VCF file (compressed or uncompressed) |
| `--rm`        | str  | Yes      | Path to the RepeatMasker `.out` file |
| `--trf`       | str  | Yes      | Path to the TRF `.dat` file |
| `-i`, `--info`| str  | Yes      | Path to the SV info file |
| `--str`       | str  | Yes      | Path to the STRchive BED file |
| `-o`, `--out` | str  | Yes      | Path to the output directory |

**Optional Arguments**

| Argument                 | Type   | Default | Description |
|--------------------------|--------|---------|-------------|
| `--min_sv_coverage`      | float  | 0.05    | Minimum intersection between a repeat element and a structural variant (0 < value < 1) |
| `--min_class_sv_coverage`| float  | 0.25    | Minimum class-level coverage to consider SV repetitive |
| `--min_total_sv_coverage`| float  | 0.75    | Minimum total repeat coverage for an SV to be considered repetitive |
| `--max_trf_overlap`      | float  | 0.1     | Maximum allowed TRF element overlap (0 < value < 1) |
| `--div`                  | float  | 0.05    | Divisor used to prioritize period size over intersection (0 < value < 1) |
| `-l`, `--len`            | int    | 100     | Diagram length (must be > 0) |
| `--rm_chain`             | str    | None    | Chain table from `rm_align_chain.py`; adds the `RM chain` block to the diagrams |
| `--debug`                | flag   | False   | Enable debug mode |
| `-h`, `--help`           | flag   | -       | Show help message and exit |


## rm_align_chain.py

Chains RepeatMasker's alignment blocks (`.align`, from RepeatMasker `-a`) into one alignment per element and writes `rm_chain.tab`, which `repeat_annotation.py --rm_chain` draws as the `RM chain` diagram block. `run_workflow.sh` runs it straight after RepeatMasker; the `.align` files are removed with the other temporary files. The table and the block are described in [Repeat/SV Diagram](Repeat-SV_diagram.md#rm-chain-block).

| Argument          | Type  | Required | Default | Description |
|-------------------|-------|----------|---------|-------------|
| `--rm_dir`        | str   | Yes      |         | Directory with the RepeatMasker `*.fa.out` and `*.fa.align` files |
| `-i`, `--info`    | str   | Yes      |         | Path to the SV info file |
| `-o`, `--out`     | str   | Yes      |         | Path to write the chain table |
| `--jump_cost`     | float | No       | 30      | Cost, in matched bases, of each jump between alignment blocks |
| `-t`, `--threads` | int   | No       | 1       | Number of RepeatMasker outputs to chain in parallel |
| `-h`, `--help`    | flag  | No       |         | Show help message and exit |

An `.out` file with hits but no `.align` is an error (RepeatMasker was run without `-a`).


## generate_plot.py

| Argument      | Required | Description                 |
|---------------|----------|-----------------------------|
| `--tsv`       | Yes      | Path to input TSV file      |
| `--out`       | Yes      | Output directory            |

## simulate_sv.py

**Required Arguments**

| Argument        | Type   | Description                                 | Default |
|-----------------|--------|---------------------------------------------|---------|
| `-m`, `--mob`   | str    | Path to the Mobile Elements file            | —       |
| `-r`, `--rep`   | str    | Path to the Repeats file                    | —       |
| `-o`, `--out`   | str    | Path to the output directory                | —       |

**Optional Arguments**

| Argument           | Type          | Description                                                                                   | Default     |
|--------------------|---------------|-----------------------------------------------------------------------------------------------|-------------|
| `--seed`           | positive int  | The seed for random generator                                                                 | 42          |
| `--len`            | positive int  | The length of the simulated base reference (before editing with SVs)                         | 100000000   |
| `-n`               | positive int  | Number of SVs to simulate                                                                     | 100         |
| `--min`            | positive int  | Minimum length of SV                                                                          | 50          |
| `--max`            | positive int  | Maximum length of SV. Not applicable to TRE and TRC. See `--read_len`                        | 50000       |
| `--flen`           | positive int  | Max detectable period size supported by TRF to determine the length of flanking sequences    | 2000        |
| `--read_len`       | positive int  | Useful to simulate TRE and TRCs that can be correctly read mapped                             | 4000        |
| `--svtypes`        | str           | File with SV types to simulate. If not provided, all types are used. See `docs/SV_simulation.md` | ""      |
| `--frac`           | flag          | Simulate SVs with fractional lengths (0.25, 0.5, 0.75) of the mobile element                  | False       |
| `--simple`         | flag          | Random pick without replacement of mobile elements and repeats                                | False       |
| `--split`          | positive int  | Distribute the BED records to n files                                                         | 1           |
| `--debug`          | flag          | Debug mode                                                                                    | False       |
| `-h`, `--help`     | flag          | Show help message and exit     

## check_vcf_ids.py

Checks, and optionally repairs, the ID column of an SV VCF. Run by `run_workflow.sh`; see
[VCF ID column](#vcf-id-column) for how the workflow uses it.

| Argument  | Required | Description |
|-----------|----------|-------------|
| `--vcf`   | Yes      | Input VCF (plain or gzip/bgzip) |
| `--mode`  | No       | `check` (default) reports missing IDs, duplicated IDs, IDs containing whitespace and records with `MATEID`, and never writes the VCF. `rewrite` writes a copy to `--out` with only the ID column changed. |
| `--out`   | `rewrite` | Output VCF. A `.gz` name is bgzip-compressed. |
| `--map`   | No       | `rewrite`: write a table of old ID to new ID here. |

Exit codes: `0` ID column usable (or nothing to rewrite), `1` not usable but repairable,
`2` not usable and not repairable, `3` input could not be read or written.

## check_ref.py

Checks that a reference FASTA is the one an SV VCF was called against. Run by
`run_workflow.sh`; see [Reference check](#reference-check) for how the workflow uses it.

| Argument   | Required | Description |
|------------|----------|-------------|
| `--vcf`    | Yes      | Input VCF (plain or gzip/bgzip) |
| `--ref`    | Yes      | Reference FASTA (plain or bgzip; indexed with `samtools faidx`, or indexable) |
| `--sample` | No       | Number of REF alleles to check (default: `200`) |

The REF check uses multi-base REF alleles (e.g. sequence-resolved DELs) when the VCF has at
least 20, and adds single-base REF alleles otherwise. An `N` on either side never counts as
a mismatch, and REF alleles that are all `N` are skipped.

Exit codes: `0` the reference matches (or nothing could be checked), `1` the reference does
not match, `3` input could not be read.
