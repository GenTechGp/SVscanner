## Creating NCI Gadi Module

The commands are just a guide only. Please do read. Please don't copy, paste and execute.

1. Create the initiating dir and file

````

DIR=/g/data/project/install

cd ${DIR}/tools
mkdir -p svclass/1.0
touch svclass/1.0/svclass

````

2. Paste the following to `svclass` file
````

#!/bin/bash

die() { echo -e "$1" >&2 ; echo ; exit 1 ; } # terminate script

module load parallel || die "could not load parallel module"

cd "SVtoolkit dir"
source "svtools/bin/activate" || die "could not activate svtools venv"
./scripts/run_classifier.sh "$@" || die "could not run_classifier.sh"

````

3. Create a module dir and module file
````
cd ${DIR}/modules
mkdir svclass
touch svclass/1.0

````

4. Paste the following to `1.0` module file
````

#%Module1.0
proc ModulesHelp { } {
    puts stderr "This module loads the svclassifier script."
}
module-whatis "SVclassifier nci bash script"

set scriptdir ${DIR}/tools/svclass/1.0
prepend-path PATH $scriptdir

````

5. Load and run module to test
````

module avail svclass
module load svclass/1.0
svclass --help

````

6. Add module path to `.bashrc` or export
````

export MODULEPATH=$DIR/modules:$MODULEPATH

````

7. Create a qsub script like below (`svclass_script.sh`). Please change parameters/paths appropriately.
````

#!/bin/bash
#PBS -P project
#PBS -N test
#PBS -l storage=gdata/kr68+gdata/ox63+gdata/if89
#PBS -l ncpus=48
#PBS -l mem=128GB
#PBS -l walltime=20:00:00
#PBS -l wd
#PBS -V

die() { echo -e "$1" >&2 ; echo ; exit 1 ; } # terminate script

export MODULEPATH=/g/data/ox63/install/modules:$MODULEPATH
module load svclass/1.0 || die "could not load svclass module"

echo "--output_dir $OUTPUT_DIR --sample $SAMPLE --sv_vcf $SV_VCF --ref_fasta $REF"
svclass --output_dir $OUTPUT_DIR --sample $SAMPLE --sv_vcf $SV_VCF --ref_fasta $REF || die "could not run svclass"

````

8. Submit the job.
````

qsub -N test -v OUTPUT_DIR=out_dir,SAMPLE=sample,SV_VCF=vcf_path,REF=ref_path svclass_script.sh

````
