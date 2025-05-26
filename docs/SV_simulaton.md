
# Structural Variant Simulation Pipeline

The following figure summarises a SV simulation pipeline comprising SV generation, embedding into a reference genome, read simulation, and variant calling.

![Illustration](/images/sv_simulation_pipeline.png)

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
| **Repeat**   | Tandem Repeat Expansion          | Repeat sequence added to base ref; BED record instructs to expand it                           |
|              | Tandem Repeat Contraction        | Repeat sequence added to base ref; BED record instructs to contract (partially delete) it      |
|              | Perfect Tandem Repetition        | No sequence added to base ref; BED record instructs VISOR to insert a perfect tandem repeat            |
|              | Approximate Tandem Repetition    | No sequence added to base ref; BED record instructs to insert an approximate tandem repeat        |

### Note
* the sequence (either Mobile or Repeat) is randomly picked and either added to the base reference or written to the BED file depending on the SV type. Then VISOR will create the SV included reference using the base reference and the BED file.

* The Mobile elements file is just a `.fasta` file (e.g. `test/databases/dfam-fasta-download.fasta`).
* The Repeat information file is expected to have four columns as given in the example file: `test/databases/GRCh38.microsatellites.bed`. However, the current implementation only requires the information in the fourth column that has the `countxmotif`.

* The argument `--svtypes` in  `scripts/sim_sv.sh` takes a file that has at most two lines (e.g., `test/databases/svtypes.txt`). The order of the two lines does not matter. If the user wants only `deletions`, `inversions` and `tandem repeat exapansions` then the file content should be like below,
```
mobile:deletion,inversion
repeat:tandem repeat expansion
```
