
# Structural Variant Simulation Pipeline

## Basic simulation
The following figure summarises a SV simulation pipeline comprising SV generation, embedding into a reference genome, read simulation, and variant calling.

![Illustration](/images/sv_simulation_pipeline.png)

## hHplotype-specific variant simulation
The following figure summarises a SV simulation pipeline specifically for diploid variants.

![Illustration](/images/sv_simulation_pipeline_diploid.png)

The bash script `scripts/sim_sv.sh` contains the pipeline.
The python script `src/simulate_sv.py` creates the [VISOR compatible SV details in bed format](https://davidebolo1993.github.io/visordoc/usage/usage.html#visor-hack).

## Structural Variant (SV) Simulation Summary

| **Sequnce Type** | **SV Type**                      | **Simulation Method**                                                                    |
|--------------|----------------------------------|------------------------------------------------------------------------------------------|
| **Mobile**   | Deletion                         | Sequence appended (added) to base ref; BED record instructs VISOR to delete it          |
|              | Inversion                        | Sequence added to base ref; BED record instructs VISOR to invert it                     |
|              | Tandem Duplication               | Sequence added to base ref; BED record instructs to make a tandem duplication           |
|              | Inverted Tandem Duplication      | Sequence added to base ref; BED record instructs to make an inverted tandem duplication |
|              | Insertion                        | No sequence added to base ref; BED record instructs VISOR to insert it                               |
| **Repeat**   | Tandem Repeat Expansion (TRE)          | Repeat sequence added to base ref; BED record instructs to expand it                           |
|              | Tandem Repeat Contraction (TRC)     | Repeat sequence added to base ref; BED record instructs to contract (partially delete) it      |
|              | Perfect Tandem Repetition (PTR)       | No sequence added to base ref; BED record instructs VISOR to insert a perfect tandem repeat            |
|              | Approximate Tandem Repetition (ATR)   | No sequence added to base ref; BED record instructs to insert an approximate tandem repeat        |

### Note
* the sequence (either Mobile or Repeat) is randomly picked and either added to the base reference or written to the BED file depending on the SV type. Then VISOR will create the SV included reference using the base reference and the BED file.

* The Mobile elements file is just a `.fasta` file (e.g. `test/databases/dfam-fasta-download.fasta`).
* The Repeat information file is expected to have four columns as given in the example file: `test/databases/GRCh38.microsatellites.bed`. However, the current implementation only requires the information in the fourth column that has the `countxmotif`.

* The argument `--svtypes` in  `src/simulate_sv.py` takes a file that has at most two lines (e.g., `test/databases/svtypes.txt`). The order of the two lines does not matter. If the user wants only `deletions`, `inversions` and `tandem repeat exapansions` then the file content should be like below,
```
mobile:deletion,inversion
repeat:tandem repeat expansion
```

* The argument `--read_len` is useful to make sure the "Repeat" type structural variants do not span longer than the average simulated read length (for TRE/PTR/ATR). This makes sure that the number of repeats in the base reference is also within the average read length (for TRC). 

* The argument `--frac` is useful to simulate SVs with fractional lengths (0.25, 0.5, 0.75) of the mobile element for deletions and insertions only.

* The argument `--split` is useful to simulate polyploid SVs by distributing svs to two or more bed files (randomly or rotation).
