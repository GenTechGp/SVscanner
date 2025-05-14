
# Structural Variant Simulation Pipeline

This following figure summarizes a simulation-based pipeline comprising SV generation, embedding into a reference genome, read simulation, and variant calling.

![Illustration](/images/sv_simulation_pipeline.png)

bash script `scripts/sim_sv.sh` contains the pipeline.
python script `src/simulate_sv.py` creates the VISOR compatible SV details in bed format.
