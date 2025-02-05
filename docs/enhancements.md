## Checker
### Limitations
#### Naming Principles 
For inversions, creating suitable names that accurately reflect the alignment pattern is challenging as they can  be misinterpreted. For inversions, the current naming convention of PASSED_LEFT and PASSED_RIGHT reflects the alignment breakpoints e.g. PASSED_RIGHT refers to the end of the alignments supporting the SV breakpoints. However, if misunderstood as specifying which SV breakpoints the read supports as it appears in the sample, the opposite definition may apply. If renamed as the latter, this may cause confusion with the REJECTED read definitions. This ambiguity is slightly improved in rejected reads by specifying _BP. Further consideration of naming may be beneficial for clarity. 

#### Primary and Supplementary alignments
In `extractReads.py`, if the current alignment is_supplementary (i.e. FLAG != 0 || FLAG != 16) the primary alignment is accessed via the SA tag. From analysing Sniffles results, the primary alignment appears to be the read that is the first entry in the string e.g. `chr19,58568863,-,11S38592M2D4277S,60,1679` in `SA:Z:chr19,58568863,-,11S38592M2D4277S,60,1679;chrX,156028611,-,37449S2190M848I2393S,60,1082;chr20,64286700,-,41293S911M23D676S,13,95;`.
 However this may not be the case for other callers, and changes to the code may be needed. But if the check for primary chr is removed for inversions and it does not significantly affect the discordance then potentially the extraction for the primary chr from the BAM file can be removed. 

#### Filtering 
* Evaluate SV caller output for additional filtering to reduce discordance 

#### Insertions/Deletions
* Breakpoint validation will not be a suitable method for INS/DEL as they do not contain the detailed read information

### Extensions
#### Insertions/Deletions
Check GT (0/1, 1/0, 1/1)


## Classifier
### Limitations/Improvements
* Need improvements to efficiency - RepeatMasker processing time is a large bottleneck (parallel)
* Simplify differences between the server tools are being run 
```
e.g. depending whether RepeatMasker is a module vs installed
find "$splitDir" -name "${sample}.*.fa" | parallel -j "$MAX_JOBS" "RepeatMasker {} -pa $THREADS_PER_JOB -html -gff -dir '$splitDir'" #NCI
find "$splitDir" -name "${sample}.*.fa" | parallel -j "$MAX_JOBS" "${repeatMasker} {} -pa $THREADS_PER_JOB -html -gff -dir '$splitDir'" #Brenner
```

ID Generation 
* Removal of custom id - test deletion of in `repeatAnnotation` `record.id=sv_id # line 524`

**Annotation of SVs** \
Ensure final annotated SVs contain all original (filtered by size) SVs even though they may not intersect a TRF or RepeatMasker entry 
* Currently we loop through all the SVs in `output_annotations()` (598) that were added based on them having an intersection with trf and rm in `read_trf()/read_rm()`. Then we search for the SV in the original VCF (520) which maybe eliminates the SVs that TRF and RM did not find any repeats for (needs double checking)
```
631 for sv_id in sv_info:
634    classification, transposition = create_vcf_record(vcf_out, sv_vcf, sv_info, sv_id, min_repetitive, strchive)

435 def create_vcf_record(vcf_file, sv_vcf, sv_info, sv_id, min_repetitive, strchive):
520     records = sv_vcf.fetch(chrom, pos-1, pos+1)
```
* Instead this could be refactored for the reverse, where you open the original SV, parse the entries and locate the annotations for that SV in the `sv_info` dict. This would also reduce additional searches by pos. 

**NON_REPETITIVE and NA labels**\
Currently 
* NON_REPETITIVE = intersection < threshold
* NA = no intersection with TRF entry or RepeatMasker Entry

Change
* NON_REPETITIVE = below threshold and not intersecting 

Most likely changes to code 
```
# Line 299-391 prepare() functions
# Change the else for classication to # 'NON_REPETITIVE'
def prepare_rm_data(rm_data):
    """Extract repeat data fields from TRF records."""
    ...
    return {
        'RM_CLASSIFICATION': ','.join(classifications) if classifications else 'NA',    
        'RM_ELEMENTS_COVERAGE': ','.join(element_coverage) if element_coverage else 'NA',
        'RM_ELEMENT_PROPORTION': ','.join(parts) if parts else 'NA',
        'RM_SV_COVERAGE': ','.join(sv_coverage) if sv_coverage else 'NA',
    }

def prepare_trf_data(trf_data):
    ...
    return {
        'TRF_CLASSIFICATION': ','.join(classifications) if classifications else 'NA',
        'TRF_PERIOD_SIZE': ','.join(period_sizes) if period_sizes else 'NA',
        'TRF_COPY_NUMBER': ','.join(copy_number) if copy_number else 'NA',
        'CONSENSUS_REPEAT': ','.join(consensus_repeats) if consensus_repeats else 'NA',
        'TRF_SV_COVERAGE': ','.join(sv_coverage) if sv_coverage else 'NA',
    }

# Also check line 468-480
if total_rm_coverage != 'NA' and float(total_rm_coverage) < min_repetitive:
    rm_data.update({key: 'NA' for key in rm_data if key != 'RM_TOTAL_SV_COVERAGE'})
    rm_data['RM_CLASSIFICATION'] = 'NON_REPETITIVE'
if total_trf_coverage != 'NA' and float(total_trf_coverage) < min_repetitive:
    trf_data.update({key: 'NA' for key in trf_data if key != 'TRF_TOTAL_SV_COVERAGE'})
    trf_data['TRF_CLASSIFICATION'] = 'NON_REPETITIVE'

```

### Extensions
#### Annotations
* Integrate phasing into pipeline


#### Visualisation 
* BED files that mark positions of repeat elements for IGV visualisation 
* Convert simple txt diagrams to interactive platform 
