# Table of Contents
1. [Installation](#installation)
2. [SVClassifier](#svclassifier)
   - [Setting up TRF + Sequence Search Engine + TE Database + RepeatMasker + Other](#setting-up-trf--sequence-search-engine--te-database--repeatmasker--other)
   - [Quick example run](#quick-example-run)
   - [Usage](#usage)
   - [Inputs & Outputs (Extended)](#inputs--outputs-extended)


## Installation 

1. Clone the repository

```
git clone https://github.com/KCCGGenomeTechLab/SVtoolkit.git
```

2. Set up Virtual Environment and install required packages. Tested with `python 3.8` and should work with higher versions as well.

```
cd SVtoolkit
python3 -m venv svtools
source svtools/bin/activate 
pip install --upgrade pip
pip install -r requirements.txt
```

## SVClassifier

### Setting up TRF + Sequence Search Engine + TE Database + RepeatMasker + Other

1. TRF (Tandem Repeats Finder), bcftools, bgzip, tabix, samtools 
```
./scripts/install_tools.sh

TRF_BINARY=$(realpath trf409.linux64)
BCFTOOLS=$(realpath bcftools-1.21/bcftools)
BGZIP=$(realpath htslib-1.21/bgzip)
TABIX=$(realpath htslib-1.21/tabix)

```

2. Download Sequence Search Engine (e.g. HMMER)
```
wget http://eddylab.org/software/hmmer/hmmer-3.4.tar.gz
tar -xvf hmmer-3.4.tar.gz
cd hmmer-3.4 
./configure
make -j 8
cd easel
make -j 8
cd ../..
rm hmmer-3.4.tar.gz

HMMER_DIR=$(realpath hmmer-3.4/src)

```
3. Download DFAM (8.9 GB - takes about 1 hour to download) and extract (67 GB - takes about 15 mins to extract)
   
   Skip if you already have downloaded a copy of the Dfam database.
   
   Then just set `DFAM_FILE` variable to the correct full path of the database file.
```
wget https://www.dfam.org/releases/Dfam_3.8/families/FamDB/dfam38-1_full.0.h5.gz
gunzip dfam38-1_full.0.h5.gz

DFAM_FILE=$(realpath dfam38-1_full.0.h5)
```
4. Download RepeatMasker
   
   Discard the error `Can't open DateRepeats: No such file or directory.` in the configuration step below.

   Make sure to run RepeatMasker on `test/sample.fa` to build libs and indices that are required for future workflows.
```
wget https://www.repeatmasker.org/RepeatMasker/RepeatMasker-4.1.6.tar.gz
tar -xvzf RepeatMasker-4.1.6.tar.gz
cd RepeatMasker
mv $DFAM_FILE Libraries/famdb
./configure -trf_prgm $TRF_BINARY -hmmer_dir $HMMER_DIR -default_search_engine hmmer
cd ..
rm RepeatMasker-4.1.6.tar.gz

REPEAT_MASKER=$(realpath RepeatMasker/RepeatMasker)
$REPEAT_MASKER test/sample.fa
```
5. Other
   
   Check if you have the following tools installed on the computer system.

   The tools can be activated on some HPCs, e.g. `module load bcftools`
```
realpath
split
parallel
```

### Quick example run

The following example processes the `test/HG002_subset_mini.vcf.gz` dataset (100 records).

It will take about 10 minutes. The majority of time is taken by the RepeatMasker step.

You may need to change the `REF_FASTA` variable to set the reference genome path (see below - Input Files).

```
./scripts/run_classifier.sh
```

### Usage

**1. Update the following variables if you are using a different dataset**


```
# Directories (change)
OUTPUT_DIR=$(realpath "SVtoolkit_output_1")

# Input Files (change)
REF_FASTA=$(realpath "/genome/hg38noAlt.fa")
SV_VCF=$(realpath "test/HG002_subset_mini/HG002_subset_mini.vcf.gz")
STR_BED=$(realpath "test/STRchive-disease-loci.bed")

# Repeat Masker and TRF programs (change if necessary)
TRF_BINARY=$(realpath "trf409.linux64")
REPEAT_MASKER=$(realpath "RepeatMasker/RepeatMasker")

BCFTOOLS=$(realpath bcftools-1.21/bcftools)
BGZIP=$(realpath htslib-1.21/bgzip)
TABIX=$(realpath htslib-1.21/tabix)

# Parameters (change if necessary)
SAMPLE="HG002_subset_mini"
NUM_SPLIT=100        
MIN_INTERSECT=0.05   
MIN_COVERAGE=0.5     
INTERVAL=0.05
DIAGRAM_LEN=100
```
**2. Adjust parameters (if necessary)**

| Parameter  | Description  | Default |
| :---- | :---- | :---- |
| `$SAMPLE` | Name of sample being analysed  |  |
| `$NUM_SPLIT` | Number of sequences/SVs in each split file  | 100 |
| `$MIN_INTERSECT` | Minimum intersection between repeat and SV to be considered before filtering | 0.05 (5%) |
| `$MIN_COVERAGE` | Minimum coverage of SV by repeats to be considered repetitive  | 0.5 (50%) |
| `$INTERVAL` | Chosen intervals to prioritise period size over intersection (for tandem repeats) | 0.05 (5%) |
| `$DIAGRAM_LEN` | Number of characters used for the diagrams | 100 |

**3. Run script**  
```
./SVclassifier/run_classifier.sh
```

### Inputs & Outputs (Extended)

#### Inputs
* `$REF_FASTA`:path to .fa file containing reference 
* `$SV_VCF`:path to .vcf file containing structural variants called by SV caller  
* `$STR_BED`:path to strchive bed file containing STR loci downloaded at https://strchive.org/loci/

#### Outputs
```
$outputDir/  
├── $sample/  
│   ├── ${sample}_annotated.vcf.gz  
│   ├── ${sample}_annotated.vcf.gz.tbi  
|   ├── ${sample}_diagrams.txt
```      
**Annotated VCF (_annotated.vcf.gz)**  
This file is the sorted and indexed VCF containing the additional annotations.  

**SV Diagrams**   
This file contains a visual representation of the repeats within the SV. The annotation to the right of the repeat under ‘SV & Flanking’ identifies the intersection. The annotation to the right of the repeat under ‘SV Diagram’ is the name of the matching repeat (RepeatMasker) or consensus repeat (TRF).

```
SV_ID chr:start_flanking-end_flanking chr1:start-end  SV_length
 SV & flanking      -------------------------[#]----------------------
 -RM-
  Simple_repeat                          	|                           	    7.5%
  Retroposon                             	|                           	    24.06%
  SINE                                   	[.]                         	    62.19%
 -TRF-
  1                                          ]                                  7.5%
 SV Diagram         [################################################]
 -RM-
  Simple_repeat   	[...]                                                       (T)n
  Retroposon          	[...........]                                  	        SVA_B
  SINE                                [..............................]       	AluYa5
 -TRF-
  1               	[...]                                       	    	    T
  ```
