## Checker
### Limitations
#### Naming Principles 
For inversions, creating suitable names that accurately reflect the alignment pattern is challenging as they can  be misinterpreted. For inversions, the current naming convention of PASSED_LEFT and PASSED_RIGHT reflects the alignment breakpoints e.g. PASSED_RIGHT refers to the end of the alignments supporting the SV breakpoints. However, if misunderstood as specifying which SV breakpoints the read supports as it appears in the sample, the opposite definition may apply. If renamed as the latter, this may cause confusion with the REJECTED read definitions. This ambiguity is slightly improved in rejected reads by specifying _BP. Further consideration of naming may be beneficial for clarity. 

#### Filtering 
* Evaluate SV caller output for additional filtering to reduce discordance 

#### Insertions/Deletions
Breakpoint validation will not be a suitable method for INS/DEL as they do not contain the detailed read information

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
* Improved id generation in extractSVs.py - currently it increments index (New ids are created incase vcfs are merged and caller recreates ids) 
* Error checking and handling for files 

### Extensions
#### Annotations
* Integrate phasing into pipeline


#### Visualisation 
* BED files that mark positions of repeat elements for IGV visualisation 
* Convert simple txt diagrams to interactive platform 
