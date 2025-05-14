
# Structural Variant Simulation Pipeline

The following figure summarises a SV simulation pipeline comprising SV generation, embedding into a reference genome, read simulation, and variant calling.

![Illustration](/images/sv_simulation_pipeline.png)

The bash script `scripts/sim_sv.sh` contains the pipeline.
The python script `src/simulate_sv.py` creates the [VISOR compatible SV details in bed format](https://davidebolo1993.github.io/visordoc/usage/usage.html#visor-hack).
