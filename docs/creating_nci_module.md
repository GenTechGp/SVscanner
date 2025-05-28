## Creating NCI Gadi Module

The commands are just a guide only. Please do read. Please don't copy, paste and execute.

1. create the initiating dir and file

````

DIR=/g/data/project/install

cd ${DIR}/tools
mkdir -p svclass/1.0
touch svclass/1.0/svclass

````

2. paste the following to `svclass` file
````

#!/bin/bash

die() { echo -e "$1" >&2 ; echo ; exit 1 ; } # terminate script

module load parallel || die "could not load parallel module"

cd "SVtoolkit dir"
source "svtools/bin/activate" || die "could not activate svtools venv"
./scripts/run_classifier.sh "$@" || die "could not run_classifier.sh"

````

3. create a module dir and module file
````
cd ${DIR}/modules
mkdir svclass
touch svclass/1.0

````

4. paste the following to `1.0` module file
````

#%Module1.0
proc ModulesHelp { } {
    puts stderr "This module loads the svclassifier script."
}
module-whatis "SVclassifier nci bash script"

set scriptdir ${DIR}/tools/svclass/1.0
prepend-path PATH $scriptdir

````

5. load and run module
````

module avail svclass
module load svclass/1.0
svclass --help

````

6. add module path to `.bashrc` or export
````

export MODULEPATH=$DIR/modules:$MODULEPATH

````