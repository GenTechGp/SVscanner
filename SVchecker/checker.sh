#!/bin/bash
#$ -S /bin/bash
#$ -terse
#$ -cwd
#$ -l mem_requested=8G 
#$ -l tmp_requested=60G 
#$ -pe smp 16
#$ -N TH126TU_checker
#$ -V


################################################################################

# Modules
export PATH=/share/ClusterShare/software/contrib/hirsam/Clair3:$PATH
export PYTHONPATH=/share/ClusterShare/software/contrib/hirsam/miniconda3/bin:$PYTHONPATH
export PATH=/directflow/KCCGGenometechTemp/projects/andmar/human_genome_pipeline/software:$PATH
module load centos7.8/phuluu/bedtools/2.29.2
# Path to virtual environment
source /directflow/KCCGGenometechTemp/projects/jasyip/checkerEnv/bin/activate

sample=sampleName
caller=Sniffles # Caller options ['Sniffles', 'CuteSV', 'SVIM']
buffer=20       # Buffer surrounding each SV breakpoint for intersection 

# Directories 
checkerDir=/directflow/KCCGGenometechTemp/projects/jasyip/SVchecker
outputDir=/directflow/KCCGGenometechTemp/projects/jasyip/analysis

# Input Files
refFASTA=/directflow/KCCGGenometechTemp/projects/andmar/human_genome_pipeline/References/hg38/hg38noAlt.fa
readsBAM=/directflow/KCCGGenometechTemp/projects/jasyip/HG002_prac/minimap/HG002.sorted.bam
svVCF=/directflow/KCCGGenometechTemp/projects/jasyip/HG002_prac/sniffles/HG002noMin.sorted.vcf.gz

################################################################################

# Output Paths
outputSampleDir=${outputDir}/${sample}
mkdir -p ${outputSampleDir}
mkdir -p ${outputSampleDir}/invResults
mkdir -p ${outputSampleDir}/dupResults
chmod +r ${outputSampleDir}/
chmod +w ${outputSampleDir}/

# Files 
readsBED=${outputSampleDir}/${sample}.bed

# DUPLICATION
# Intermediate Files
dupBED=${outputSampleDir}/dupResults/dup_caller.bed
dupSupport=${outputSampleDir}/dupResults/dup_caller_support.tab
dupBuffer=${outputSampleDir}/dupResults/dup_caller_buffer.bed
dupDepth=${outputSampleDir}/dupResults/dup_depth.tab
dupOverlaps=${outputSampleDir}/dupResults/dup_overlaps.tab
dupOrdered=${outputSampleDir}/dupResults/dup_ordered.tab
dupSupportingReads=${outputSampleDir}/dupResults/dup_supporting_reads.tab
# dupBED=$(mktemp)
# dupSupport=$(mktemp)
# dupBuffer=$(mktemp)
# dupDepth=$(mktemp)
# dupOverlaps=$(mktemp)
# dupOrdered=$(mktemp)
# dupSupportingReads=$(mktemp)

# Final Outputs
dupSupportingReadsDetails=${outputSampleDir}/dupResults/dup_supporting_read_details.tab
dupDiscordant=${outputSampleDir}/dupResults/dup_discordant.tab
dupChecked=${outputSampleDir}/dupResults/dup_checked.tab

# INVERSION 
# Intermediate Files
invBED=${outputSampleDir}/invResults/inv_caller.bed
invSupport=${outputSampleDir}/invResults/inv_caller_support.tab
invBuffer=${outputSampleDir}/invResults/inv_caller_buffer.bed
invDepth=${outputSampleDir}/invResults/inv_depth.tab
invOverlaps=${outputSampleDir}/invResults/inv_overlaps.tab
invOrdered=${outputSampleDir}/invResults/inv_ordered.tab
invSupportingReads=${outputSampleDir}/invResults/inv_supporting_reads.tab
# invBED=$(mktemp)
# invSupport=$(mktemp)
# invBuffer=$(mktemp)
# invDepth=$(mktemp)
# invOverlaps=$(mktemp)
# invOrdered=$(mktemp)
# invSupportingReads=$(mktemp)

# Final Outputs 
invSupportingReadsDetails=${outputSampleDir}/invResults/inv_supporting_read_details.tab
invDiscordant=${outputSampleDir}/invResults/inv_discordant.tab
invChecked=${outputSampleDir}/invResults/inv_checked.tab

## CUSTOM PROGRAMS
BAM2BED=${checkerDir}/extractReads.py
MERGE=${checkerDir}/mergeTables.py
DISCORDANT=${checkerDir}/discordant.py
INV_DEPTH=${checkerDir}/invDepth.py
DUP_DEPTH=/${checkerDir}/dupDepth.py
ORDER=${checkerDir}/orderAlignments.py
INVERSION=${checkerDir}/invCheck.py
DUPLICATION=${checkerDir}/dupCheck.py

######################## Caller Variable Definitions ###########################
if [ "$caller" = "Sniffles" ]; then # sniffles -h
    mapq=25
    min_alignment_length=1000
    max_splits_base=3
    max_splits_kb=0.1
    svQuery='%CHROM\t%POS\t%INFO/END\t%ID\t%INFO/SVTYPE\t%INFO/RNAMES\t%FILTER\n'
elif [ "$caller" = "cuteSV" ]; then # https://github.com/tjiangHIT/cuteSV 
    mapq=20
    min_alignment_length=500
    max_splits_base=7
    max_splits_kb=0 
    svQuery='%CHROM\t%POS\t%INFO/END\t%ID\t%INFO/SVTYPE\t%INFO/RNAMES\t%FILTER\n'
elif [ "$caller" = "SVIM" ]; then # svim reads --help
    mapq=20
    # Note: SVIM does not filter based on the below (low thresholds are set)
    min_alignment_length=100
    max_splits_base=7
    max_splits_kb=0.1 
    svQuery='%CHROM\t%POS\t%INFO/END\t%ID\t%INFO/SVTYPE\t%INFO/RNAMES\t\n'
else # Use default values of Sniffles
    mapq=25
    min_alignment_length=1000
    max_splits_base=3
    max_splits_kb=0.1 
fi

########################## Alignment Preprocessing #############################
# 1. Convert bam using external program with the FLAGS and alignment info
python $BAM2BED -i $readsBAM > $readsBED
# chr   start   end    readID   MAPQ    flag    mismatch primary_tag    SA_tag

####################### Structural Variant Preprocessing #######################

## 1. Get breakpoints of Sniffles (VCF) output for duplications in single bed format 
# Extracts location, id, sv and supporting read information and outputs to files 
process_svtype() {
    local svType="$1"       # SVTYPE (e.g., DUP or INV)
    local svBED="$2"        # Output BED file
    local svSupport="$3"    # Output support file

    # Process the VCF with bcftools and awk
    bcftools query -f "${svQuery}" "${svVCF}" | \
    awk -v svtype="${svType}" -v svBED="${svBED}" -v svSupport="${svSupport}" 'OFS="\t" {
        if ($1 ~ /^chr([1-9]$|1[0-9]$|2[0-2]$|X$|Y$)/ && $5 == svtype) {
            # Count the number of RNAMES by splitting them on commas
            # NOTE: Count instead of using INFO:SUPPROT/RE as the flag differs between callers 
            numSupport = split($6, rnames, ",")  

            # If $7 (FILTER-SVIM) is missing, output "NA" instead
            filterField = ($7 == "" ? "NA" : $7)

            # CHR,POS,END,ID - remove SVTYPE column ($5)
            print $1, $2, $3, $4 > svBED

            # ID, SUPPORT (count of RNAMES), RNAMES, FILTER
            print $4, numSupport, $6, $7 > svSupport
        }
    }'
}

# Process different SVTYPEs
process_svtype "DUP" "$dupBED" "$dupSupport"
process_svtype "INV" "$invBED" "$invSupport"

################################# Duplications #################################
if [[ ! -s ${dupBED} || ! -s ${dupSupport} ]]; then
    echo "File (${svVCF}) does not contain any duplications"
else
    ## 2. Transform sniffles bed into start and end bed files with buffer --> Convert so start and end positions are separated 

    awk -v buffer=${buffer} 'OFS="\t" {
        startL = $2 - buffer
        if (startL < 1) startL = 1
        endL = $2 + buffer

        startR = $3 - buffer
        if (startR < 1) startR = 1
        endR = $3 + buffer

        print $1, startL, endL, $4".L"  
        print $1, startR, endR, $4".R"
    }' ${dupBED} > ${dupBuffer}
    
    ## 3. Get the Depth Info
    python3 ${DUP_DEPTH} -bam $readsBAM -bed $dupBED -ref $refFASTA -o $dupDepth

    ## INTERSECT STRUCTURAL VARIANT WITH READS (BED)
    bedtools intersect -wa -wb -a $dupBuffer -b $readsBED | \
    awk 'OFS="\t" {print $1, $2, $3, $4, $6, $7, $8, $9, $10, $11, $12, $13, $14}' > $dupOverlaps   

    # ORDER SUPPLEMENTARY ALIGNMENTS
    python3 $ORDER -i $dupOverlaps -o $dupOrdered

    # RUN READ CHECKER based on caller default 
    if [ "$caller" = "Sniffles" ]; then 
        python3 $DUPLICATION -i $dupOrdered -d $dupDepth -o $dupSupportingReads -r $dupSupportingReadsDetails -mapq 25 -min-alignment-length 1000 -max-splits-kb 0.1 -max-splits-base 3 
    elif [ "$caller" = "cuteSV" ]; then
        python3 $DUPLICATION -i $dupOrdered -d $dupDepth -o $dupSupportingReads -r $dupSupportingReadsDetails -mapq 20 -min-alignment-length 500 -max-splits-kb 0 -max-splits-base 7 
    else
        python3 $DUPLICATION -i $dupOrdered -d $dupDepth -o $dupSupportingReads -r $dupSupportingReadsDetails

    python3 ${DISCORDANT} -c $dupSupportingReads -i $dupSupport -o $dupDiscordant

    ### DUPLICATION COMBINE
    # dupSupprotingReads
    # $1 $2             $3                 $4  $5   $6  $7  $8          $9          
    # ID #CALLER_PASSED #SupportingReads chr start end depth left_depth right_depth left_ratio  right_ratio 
    columns="ID\tchr\tstart\tend\t\
    depth\tleftDepth\trightDepth\tleftDepthRatio\trightDepthRatio\tCALLER_SUPPORT\t\
    TOTAL_PASSED\tTANDEM_DUPLICATION\tTANDEM_REPEAT\t\
    INTERSPERSED_DUPLICATION\tINTERSPERSED_REPEAT\t\
    TOTAL_REJECTED\tMISSING_INTERSECT\tMISSING_PAIRS\tFILTER"
    echo -e "$columns" > $dupChecked
    python3 ${MERGE} -c 0  $dupSupportingReads $dupSupport | awk 'OFS="\t" {print $1, $4, $5, $6, $7, $8, $9, $10, $11, $19, $2, $12, $13, $14, $15, $16, $17, $18, $21}' >> $dupChecked
fi

################################## Inversions ##################################
if [[ ! -s ${invBED} || ! -s ${invSupport} ]]; then
    echo "No <INV> to check in $svVCF"
else
    awk -v buffer=${buffer} 'OFS="\t" {
        print $1, $2-buffer, $2+buffer, $4".L"
        print $1, $3-buffer, $3+buffer, $4".R" 
    }' ${invBED} > ${invBuffer}

    ## 3. Get the Depth Info
    python3 ${INV_DEPTH} -bam $readsBAM -bed $invBED -ref $refFASTA -o $invDepth

    bedtools intersect -wa -wb -a $invBuffer -b $readsBED | \
    awk 'OFS="\t" {print $1, $2, $3, $4, $6, $7, $8, $9, $10, $11, $12, $13, $14}' > $invOverlaps  
    
    # ORDER SUPPLEMENTARY ALIGNMENTS
    python3 $ORDER -i $invOverlaps -o $invOrdered

    # RUN READ CHECKERS
    python3 $INVERSION -i $invOrdered -d $invDepth -o $invSupportingReads -r $invSupportingReadsDetails
    python3 ${DISCORDANT} -c $invSupportingReads -i $invSupport -o $invDiscordant

    ### INVERSION COMBINE
    columns="ID\tchr\tstart\tend\t\
    startDepth\tendDepth\t\CALLER_SUPPORT\t\
    TOTAL_PASSED\tEXTENSION\tMULTI\t\
    TOTAL_REJECTED\tLEFT_BP\tRIGHT_BP\tBOTH\tORIENTATION\tFILTER" 
    echo -e "$columns" > $invChecked                                                                          
    python3 ${MERGE} -c 0  $invSupportingReads $invSupport | awk 'OFS="\t" {print $1, $4, $5, $6, $7, $8, $16, $2, $9, $10, $11, $12, $13, $14, $15, $18}' >> $invChecked

    echo "Checked Inversions: $invChecked"
    echo "Discordant Reads: $invDiscordant"
    echo "Read Details: $invSupportingReadsDetails"
fi

####################### Cleanup: Remove Intermediates ##########################
# filesToRemove=(
#     "$invBED" "$invSupport" "$invBuffer" "$invDepth" "$invOverlaps" "$invOrdered" "$invSupportingReads"
#     "$dupBED" "$dupSupport" "$dupBuffer" "$dupDepth" "$dupOverlaps" "$dupOrdered" "$dupSupportingReads"
# )
# # Loop through and remove each file
# for file in "${filesToRemove[@]}"; do
#     rm -f "$file"
# done












