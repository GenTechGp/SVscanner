#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ---------------------------
// Parameters
// ---------------------------
params.vcf     = null
params.ref     = null
params.outdir  = 'results'
params.species = 'human'
params.str_bed = ''
params.cpus    = 1

// ---------------------------
// Process
// ---------------------------
process SVSCANNER_ANNOTATE {

    tag { vcf.baseName }
    cpus params.cpus

    // NO conda or self containerization
    // We rely on `module load` externally

    publishDir params.outdir, mode: 'copy', overwrite: true, pattern: 'svscanner.annotated.vcf.gz*'

    input:
      tuple path(vcf), path(ref), val(species), val(str_bed)

    output:
      path "svscanner.annotated.vcf.gz"
      path "svscanner.annotated.vcf.gz.tbi"

    script:
      def strArg = str_bed ? "--str_bed ${str_bed}" : ""

      """
      set -euo pipefail

      echo "==== SVscanner Nextflow wrapper ===="
      echo "VCF: ${vcf}"
      echo "REF: ${ref}"
      echo "Species: ${species}"
      echo "Working dir: \$(pwd)"

      # Prepare sandbox directories
      mkdir -p scripts src test/databases

      # Copy workflow scripts/resources
      cp -r "${projectDir}/../scripts/"* scripts/
      cp -r "${projectDir}/../src/"* src/

      # Default STR file
      if [ -z "${str_bed}" ]; then
        cp "${projectDir}/../test/databases/STRchive-disease-loci.bed" test/databases/STRchive-disease-loci.bed
      fi

      # Run the SVScanner workflow
      bash "scripts/run_workflow.sh" \\
        --out "svscanner_out" \\
        --vcf "${vcf}" \\
        --ref "${ref}" \\
        --species "${species}" \\
        ${strArg} \\

      # Validate and publish output
      if [ ! -f "svscanner_out/annotated.vcf.gz" ]; then
        echo "ERROR: Expected svscanner_out/annotated.vcf.gz not found" >&2
        exit 1
      fi

      mv "svscanner_out/annotated.vcf.gz" svscanner.annotated.vcf.gz

      if [ -f "svscanner_out/annotated.vcf.gz.tbi" ]; then
        mv "svscanner_out/annotated.vcf.gz.tbi" svscanner.annotated.vcf.gz.tbi
      else
        tabix -p vcf svscanner.annotated.vcf.gz
      fi

      echo "Annotation complete."
      """
}


// ---------------------------
// Workflow
// ---------------------------
workflow {

    if (!params.vcf) error "Missing required param: --vcf"
    if (!params.ref) error "Missing required param: --ref"

    input_ch = Channel.of(
      tuple(
        file(params.vcf),
        file(params.ref),
        params.species,
        params.str_bed,
      )
    )

    SVSCANNER_ANNOTATE(input_ch)
}
