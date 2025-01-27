# Table of Contents
1. [Installation](#installation)
2. [SVChecker](#svchecker)
   - [Usage](#usage)
   - [Inputs & Outputs (Extended)](#inputs--outputs-extended)
3. [SVClassifier](#svclassifier)
   - [Setting up TRF and RepeatMasker](#setting-up-trf-and-repeatmasker)
   - [Usage](#usage)
   - [Inputs & Outputs (Extended)](#inputs--outputs-extended)


## Installation 

1. Download repo

```
git clone https://github.com/KCCGGenomeTechLab/SVtoolkit.git
```

2. Set up Virtual Environment

```
python3 -m venv svtools
source svtools/bin/activate 
pip install --upgrade pip  
```

3. Install requirements
```
pip install -r requirements.txt
```

## SVChecker

### Usage

**1. Update paths to input/output files (examples given)**
```
# Load bedtools, samtools, bcftools
export PATH=...

# Specify path to activate virtual environment  
source ../svtools/bin/activate

sample=HG002  
caller=Sniffles         # Caller options ['Sniffles', 'CuteSV', 'SVIM'] 
buffer=20               # Buffer surrounding each SV breakpoint for intersection

# Directories
checkerDir=.../SVtools/SVchecker
outputDir=.../outputDir # Specify a directory for output files 

## Input Files  
refFASTA=.../hg38noAlt.fa  
readsBAM=.../HG002.sorted.bam 
svVCF=.../HG002.sorted.vcf.gz
```

**2. Adjust parameters**

| Parameter | Description  | Default |
| ----- | ----- | ----- |
| `$sample` | Name of sample being analysed  |  |
| `$buffer` | buffer added to the start (POS-buffer) and end (END+buffer) of SV breakpoints to account for technical variation | 20 |
| `$caller` | Identify caller used to generate VCF [Sniffles, CuteSV, SVIM] | Sniffles |
| Supporting Reads Filtering (declared based on `$caller`) |  |  |
| `$mapq` | Minimum mapping quality value of alignment to be taken into account | Caller dependent (e.g. 25) |
| `$min_alingment_length` | Reads with alignments shorter than this length (in bp) will be ignored  | Caller dependent (e.g. 1000) |
| `$max_splits_base` | Base number of split segments (alignments) a read may be aligned before it is ignored | Caller dependent (e.g. 3) |
| `$max_splits_kb` | Additional number of splits per kilobase read sequence allowed before reads are ignored  | Caller dependent (e.g. 0.1) |

**3. Run script**
```
chmod +x checker.sh
./checker.sh 
```

### Inputs & Outputs (Extended)

#### Input

* `$readsBAM` - path to bam file containing reads alignmed by minimap2  
* `$svVCF` - path to vcf file containing structural variants called by SV caller  
* `$refFASTA` - path to .fa file containing reference 

#### Output
```
$outputDir/  
├── $sample/  
│   ├── invResults/  
|   |   ├── inv_checked.tab  
|   |   ├── inv_discordant.tab  
|   |   ├── inv_supporting_read_details.tab  
│   ├── dupResults/  
|   |   ├── dup_checked.tab  
|   |   ├── dup_discordant.tab  
|   |   ├── dup_supporting_read_details.tab  
```  

**SV Summary (_checked.tab)**

This file includes the depth of coverage at the breakpoints (inversions) or an average of flanking regions (duplications), a comparison between the number of supporting reads from the caller and the checker, a count of reads that met flags specific to the SV type and whether the SV was filtered by the caller.

| svID | chr | start | end | depth | # Caller Support | # Checker Support | Supporting Read Flags Count | Rejected Read Flags Count | Filter |
| :---- | :---- | :---- | :---- | :---- | :---- | :---- | :---- | :---- | :---- |


**Discordant Reads (_discordant.tab)**  
This file lists the number of supporting reads by the caller and checker with the ids of discordant reads.

| svID | # Caller Support | # Checker Support | # Discordant  | Discordant RNAMES |
| :---- | :---- | :---- | :---- | :---- |
| INV.18A8AS0   | 8 | 7 | 1  | fed0560b |
| INV.C8EBS1 | 8 | 8 | 0 | NA |

**Supporting Reads Details (supporting_read_details.tab)**  
This file includes the flag of every read for each SV.

| svID | readID | flag |
| :---- | :---- | :---- |
| INV.C922S0  | 285b70ae | PASSED_RIGHT |
| INV.C922S0  | 940a1765 | MISSING_BOTH |

## SVClassifier

### Setting up TRF and RepeatMasker

#### TRF (Tandem Repeats Finder) 
```
wget https://github.com/Benson-Genomics-Lab/TRF/releases/download/v4.09.1/trf409.linux64
chmod +x trf409.linux64
TRF_PRG=$PWD/trf409.linux64
```

#### RepeatMasker 
1. Download Sequence Search Engine (e.g. HMMER)
```
wget http://eddylab.org/software/hmmer/hmmer-3.4.tar.gz
tar -xvf hmmer-3.4.tar.gz
./hmmer-3.4/configure
make -C hmmer-3.4 -j 8
make -C hmmer-3.4/easel -j 8
HMMER_DIR=$PWD/hmmer-3.4/src
```
2. Download DFAM (8.9 GB)
   
   Skip if you have already downloaded a copy of the Dfam database.
   
   Then just set `DFAM_FILE` variable to the correct full path of the database file.
```
wget https://www.dfam.org/releases/Dfam_3.8/families/FamDB/dfam38-1_full.0.h5.gz
gunzip dfam38_full.0.h5.gz
DFAM_FILE=$PWD/dfam38_full.0.h5
```
3. Download RepeatMasker
```
wget https://www.repeatmasker.org/RepeatMasker/RepeatMasker-4.1.6.tar.gz
tar -xvzf RepeatMasker-4.1.6.tar.gz
cd RepeatMasker
mv $DFAM_FILE Libraries/famdb
./configure -trf_prgm $TRF_PRG -hmmer_dir $HMMER_DIR
cd ..
```

### Usage

**1. Update paths to input/output files (examples given)**


```
# Load bedtools, samtools, bcftools
export PATH=...
source ../svtools/bin/activate # Activate virtual environment

sample=HG002  
numSplit=100        # Number of sequences per file

# Directories
classifierDir=.../SVtools/SVchecker

outputSampleDir=${outputDir}/${sample}  
splitDir=${outputSampleDir}/${sample}_${numSplit}

# Input Files  
refFASTA=.../hg38noAlt.fa
svVCF=.../HG002.sorted.vcf.gz

# Repeat Masker and TRF programs 
TRF_BINARY=.../trf409.linux64  
repeatMasker=.../RepeatMasker  
export PATH=.../hmmer-3.4/src:$PATH
```
**2. Adjust parameters**

| Parameter  | Description  | Default |
| :---- | :---- | :---- |
| `$sample` | Name of sample being analysed  |  |
| `$numSplit` | Number of sequences/SVs in each split file  | 100 |
| `$minIntersect` | Minimum intersection between repeat and SV to be considered before filtering | 0.05 (5%) |
| `$minCoverage` | Minimum coverage of SV by repeats to be considered repetitive  | 0.5 (50%) |
| `$interval` | Chosen intervals to prioritise period size over intersection (for tandem repeats) | 0.05 (5%) |
| `$diagramLength` | Number of characters used for the diagrams | 100 |

**3. Run script**  
```
chmod +x classifier.sh
./classifier.sh 
```

### Inputs & Outputs (Extended)

#### Inputs
* `$svVCF`   	path to .vcf file containing structural variants called by SV caller  
* `$refFASTA` 	path to .fa file containing reference 

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
