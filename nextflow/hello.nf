#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Default parameter
params.msg = params.msg ?: "Hello, Nextflow DSL2 World!"

// A simple process that prints the message
process HELLO {
    input:
    val message

    output:
    stdout

    script:
    """
    echo "${message}"
    """
}

// Workflow definition
workflow {
    // Pass the parameter to the process
    HELLO(params.msg)
}
