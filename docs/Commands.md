# Commands
1. [extrract_sv.py](#extrract_svpy)
2. [repeat_annotatoin.py](#repeat_annotatoinpy)
3. [generate_plot.py](#generate_plotpy)
3. [simulate_sv.py](#simulate_svpy)

## extrract_sv.py

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
| `--debug`                | flag   | False   | Enable debug mode |
| `-h`, `--help`           | flag   | -       | Show help message and exit |


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

