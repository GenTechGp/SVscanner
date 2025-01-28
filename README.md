# SV Toolkit

*A collection of tools for validating and annotating structural variants (SV).* 

## Overview

**SVChecker** is a tool that assists in validating SV calls by reviewing supporting evidence provided by SV callers. Via a script, the tool takes as input the reads (BAM), called variants (VCF) and reference genome (FASTA). The checker evaluates various types of evidence, such as reads overlapping breakpoints, read orientation, and coverage, to determine whether any given variants require further review.

**SVClassifier** is a tool for annotating tandem repeats and mobile elements identified by Tandem Repeat Finder (TRF) and RepeatMasker within SVs. Via a script, the tool takes as input the file containing called variants (VCF) and reference genome (FASTA). It runs instances of TRF and RepeatMasker, then identifies and prioritises entries intersecting with the SV for annotation or repetitive characteristics. Additionally, the tool enables visualisation of the annotation entries via a text file. 

## Workflow

#### SVChecker  
![GS](/images/checker_workflow_dark.svg#gh-dark-mode-only)
![GS](/images/checker_workflow.svg#gh-light-mode-only)

#### SVClassifier
![](/images/classifier_workflow.svg) 


### Main Requirements

#### Python Requirements 
* Python 3.8+ required  
* `pysam`
* `numpy `
* `pandas`
* `tqdm`
* `h5py` (for RepeatMasker)

#### Other
* `parallel` installed on system

#### SVClassifier

Executable versions of annotation software: 
* Tandem Repeat Finder (https://github.com/Benson-Genomics-Lab/TRF/releases/tag/v4.09.1)  
* RepeatMasker (https://www.repeatmasker.org/RepeatMasker/) or see [Usage.md](https://github.com/KCCGGenomeTechLab/SVtoolkit/blob/main/docs/usage.md#RepeatMasker)


### Main Input files 

#### SVChecker

* Aligned, sorted bam containing long read alignments from **minimap2**  
* VCF (vcf.gz) containing called SVs from **Sniffles2, cuteSV, SVIM** (must contain INFO flag `'RNAMES'`)  
  * SVIM must be run with `--read_names` (`svim reads --help`)

#### SVClassifier

* VCF (.vcf.gz) containing called SVs (must contain INFO flag `'SVTYPE'`)  
  * ***NOTE***: Does not assess symbolic insertions, translocations
* *Optional*: STRchive dataset for genome (in extended BED format)
  * [hg19](https://strchive.org/_astro/STRchive-disease-loci.hg19.DWACvaXd.bed) [hg38](https://strchive.org/_astro/STRchive-disease-loci.hg38.DR-UScgX.bed) [T2T-chm13](https://strchive.org/_astro/STRchive-disease-loci.T2T-chm13.Cm-HAugT.bed)

### Main Output files 

#### SVChecker

* **SV Summary** - TXT file containing summary of SVs and their supporting reads for inversions and duplications `(_checked.tab)`
* **Discordant Reads** - TXT file containing discordant reads between caller and checker for inversions and duplications `(..._discordant.tab)`
* **Read Details** - TXT file containing details of reads `(..._supporting_read_details.tab)`

#### SVClassifier

* VCF annotated with TRF and RepeatMasker info  
* TSV containing TRF and RepeatMasker repeats with extended detail
* TXT file containing diagrams of SV and repeats 

## How to Run 

Details on how to [run on SGE and extended details for input/output files](/docs/usage.md) 

## User Guide (Algorithms)

Explanation for [SVChecker](/docs/user_guide.md) algorithms  
Explanation for [SVClassifier](/docs/user_guide.md) algorithms
