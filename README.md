# SV Toolkit

*A workflow to annotate and analyse structural variants (SV) with Repeat information.* 

## Overview

A workflow for annotating tandem repeats and mobile elements identified by Tandem Repeat Finder (TRF) and RepeatMasker (RM) within SVs. The workflow takes SV information (VCF) and reference genome (FASTA) as inputs. First, flanking sequences around the SV are extracted. Then extracted sequences are annotated using TRF and RM. The workflow has three main outputs.

1. SV VCF file with repeat information annotated.
2. Diagrams with Repeat/SV annotations.
3. Histograms and density plots summarising Repeat/SV annotations.

![Illustration](/images/SVscanner_repeat_annotation_workflow.png)

- Read more about the SV extraction step [here](docs/extract_sv_and_seq_consensus_documentation.md)
- Read more about the repeat annotation step [here](docs/repeat_annotation_steps.md)
- Check commandline arguments for [extract_sv.py](), [repeat_annotation.py](), [generate_plots.py]()

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

* VCF (.vcf.gz) containing called SVs (must contain INFO flag `'SVTYPE'`)  
  * ***NOTE***: Does not assess symbolic insertions, translocations
* *Optional*: STRchive dataset for genome (in extended BED format)
  * [hg19](https://strchive.org/_astro/STRchive-disease-loci.hg19.DWACvaXd.bed) [hg38](https://strchive.org/_astro/STRchive-disease-loci.hg38.DR-UScgX.bed) [T2T-chm13](https://strchive.org/_astro/STRchive-disease-loci.T2T-chm13.Cm-HAugT.bed)

### Main Output files 

#### SVClassifier

* VCF annotated with TRF and RepeatMasker info  
* TSV containing TRF and RepeatMasker repeats with extended detail
* TXT file containing diagrams of SV and repeats 

## How to Run 

Details on how to [run on SGE and extended details for input/output files](/docs/usage.md) 

