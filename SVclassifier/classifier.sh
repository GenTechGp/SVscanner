#!/bin/bash
#$ -S /bin/bash
#$ -terse
#$ -cwd
#$ -l mem_requested=16G 
#$ -l tmp_requested=60G 
#$ -pe smp 32
#$ -N classifierHG002
#$ -V

################################################################################

# Modules 
export PATH=/share/ClusterShare/software/contrib/hirsam/Clair3:$PATH
export PYTHONPATH=/share/ClusterShare/software/contrib/hirsam/miniconda3/bin:$PYTHONPATH
export PATH=/directflow/KCCGGenometechTemp/projects/andmar/human_genome_pipeline/software:$PATH
module load centos7.8/phuluu/bedtools/2.29.2
# Path to virtual environemnt 
source /directflow/KCCGGenometechTemp/projects/jasyip/checkerEnv/bin/activate

sample=sampleName
numSplit=100        # Number of sequences per file 
minIntersect=0.05
minCoverage=0.5
interval=0.05
diagramLength=100

# Directories
classifierDir=/directflow/KCCGGenometechTemp/projects/jasyip/SVclassifier
outputDir=/directflow/KCCGGenometechTemp/projects/jasyip/performance
outputSampleDir=${outputDir}/${sample}
splitDir=${outputSampleDir}/${sample}_${numSplit}

# Input Files
refFASTA=/directflow/KCCGGenometechTemp/projects/andmar/human_genome_pipeline/References/hg38/hg38noAlt.fa
svVCF=/directflow/KCCGGenometechTemp/projects/jasyip/analysis/ANSH56_RF/KISKUM_ataxia.ANSH56_RF.hg38.SVs.phased.vcf.gz

# Repeat Masker and TRF programs
TRF_BINARY=/directflow/KCCGGenometechTemp/projects/iradev/tandem_repeats/scripts/trf409.linux64
repeatMasker=/directflow/KCCGGenometechTemp/projects/jasyip/software/RepeatMasker/RepeatMasker
export PATH=/directflow/KCCGGenometechTemp/projects/jasyip/software/hmmer-3.4/src:$PATH

################################################################################

# Change locale settings 
export LANGUAGE=en_US.UTF-8
export LC_ALL=en_US.UTF-8
export LANG=en_US.UTF-8

# Setup allocation for RepeatMasker parallelisation
NTHREADS=$NSLOTS    # Total number of threads
MAX_JOBS=8          # Max number of RepeatMasker process to run in parallel 
THREADS_PER_JOB=$((NTHREADS / MAX_JOBS)) # Number of threads allocated to each RepeatMasker job (internal)

# File Intermediates
svTAB=${outputSampleDir}/variant_flanking.tab
IDfile=${outputSampleDir}/${sample}_id.tab
RMfile=${outputSampleDir}/${sample}_rm.tab
TRFfile=${outputSampleDir}/${sample}_trf.tab
# svTAB=$(mktemp)
# IDfile=$(mktemp)
# RMfile=$(mktemp)
# TRFfile=$(mktemp)

# Final Outputs
visOutput=${outputSampleDir}/${sample}_diagrams.txt
annotatedVCF=${outputSampleDir}/${sample}_annotated.vcf
rmTSV=${outputSampleDir}/${sample}_annotatedRM.tsv
trfTSV=${outputSampleDir}/${sample}_annotatedTRF.tsv

# CUSTOM PROGRAMS
VARIANT_FLANKING=${classifierDir}/extractSVs.py
TO_FASTA=${classifierDir}/toFasta.sh
ANNOTATE=${classifierDir}/repeatAnnotation.py
VISUALISE=${classifierDir}/repeatDiagram.py
# PARALLEL_TO_FASTA=${classifierDir}/parallelToFasta.sh

####################### Structural Variant Preprocessing ##########################
mkdir -p $splitDir
cd ${outputSampleDir}

# 1) Extract sequence and flanking regions for variants
echo "1. Extracting structural variant sequences from VCF"
python3 $VARIANT_FLANKING -vcf $svVCF -fa $refFASTA -out $svTAB
echo ..Done

# 2) Create fasta sequences
echo "2. Converting structural variant sequences to FASTA"
split --numeric-suffixes=1 --suffix-length=3 -l $numSplit "$svTAB" $splitDir/${sample}.
chmod +r $splitDir/*
chmod +w $splitDir/*
find "$splitDir" -name "${sample}.*" | parallel -j $NTHREADS "
    output_fasta=\"{}.fa\"
    id_file=\"{}_id.tab\"
    bash \"${TO_FASTA}\" {} \"\$output_fasta\" \"\$id_file\""
echo ..Done

############################ Annotation Programs ###############################
# 3) Run Tandem Repeat Finder and RepeatMasker - wait for both to complete
echo "3. Running Tandem Repeat Finder"
find "$splitDir" -name "${sample}.*.fa" | parallel -j $NTHREADS --bar "${TRF_BINARY} {} 2 7 7 80 10 50 500 -h -ngs > {.}.dat"
echo ..Done
echo "3. Running RepeatMasker"
find "$splitDir" -name "${sample}.*.fa" | parallel -j $MAX_JOBS "${repeatMasker} {} -pa $THREADS_PER_JOB -html -gff -dir ${splitDir}"
wait 
echo "..Done"

# Process the RepeatMasker files (remove the header and 16th column)
for rm_output in "$splitDir"/*.fa.out; do
    tail -n +4 "$rm_output" | awk '{
        if ($16 == "*") $16 = "";
        else $16 = "";
        print $0
    }' OFS='\t' > "${rm_output}.tab"
done

# ############################### Classify SVs ###################################
# # COMBINE SPLIT FILES
cat ${splitDir}/*_id.tab > $IDfile
cat ${splitDir}/*.out.tab > $RMfile
cat ${splitDir}/*.dat > $TRFfile

echo "4. Running visualiser"
python3 ${VISUALISE} -sv $IDfile -trf $TRFfile -rm $RMfile -out $visOutput -length $diagramLength -min $minIntersect

echo "5. Annotating SV VCF with repeat information"
python3 ${ANNOTATE} -id $IDfile -trf $TRFfile -rm $RMfile -vcf $annotatedVCF -trf_tsv $trfTSV -rm_tsv $rmTSV -sv_vcf $svVCF -min $minIntersect -mr $minCoverage -div $interval

# Sort and Index the annotated VCF
bcftools sort -Oz -o ${annotatedVCF}.gz ${annotatedVCF}
bcftools index -t ${annotatedVCF}.gz 

####################### Cleanup: Remove Intermediates ##########################
# rm -r $splitDir
# filesToRemove=("$svTAB" "$IDfile" "$RMfile" "$TRFfile")
# # Loop through and remove each file
# for file in "${filesToRemove[@]}"; do
#     rm -f "$file"
# done


