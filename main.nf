#!/usr/bin/env nextflow
nextflow.enable.dsl=2

start_var = Channel.from("""
*********Start running Metavir pipeline*********
nf-metavir is a workflow for metagenomics qc, assembly and taxonomic classification
**************************************
""")
start_var.view()

if( !nextflow.version.matches('20.+') ) {
    ch_ver=Channel.from("This workflow requires Nextflow version 20.07 or greater -- You are running version $nextflow.version").view()
    exit 1
}

// Help Message
def helpMSG() {
    log.info """
    *********Assembly and taxonomic classification workflow for (viral) metagenomics*********

        Usage example:
    nextflow run main.nf --illumina illumina/ --assembler megahit -profile planet
    or
    nextflow run main.nf --illumina illumina/ --assembler megahit --host_map host_ref/ -profile uppmax
        Input:
    --illumina                  path to the directory containing the illumina read file (fastq) (default: $params.illumina)
    --assembler                 the assembler to use in the assembly step (default: $params.assembler)
        Optional input:
    --host_ref                  path to the host reference genome to map on
    --k2nt_db                   path to the Kraken2 nucleotide database (e.g. nt)
    --k2prot_db                 path to the Kraken2 protein database (e.g. nr)
    --diamond_db                path to the Diamond protein database (e.g. nr)

        Output:
    --output                    path to the output directory (default: $params.output)

        Outputed files:
    qc                          The reads file after qc, qc logs and host mapping logs
    reads                       The taxonomic classifications at reads level
    contigs                     The taxonomic classifications at contigs level

        Basic Parameter:
    --cpus                      max cores for local use [default: $params.cpus]
    --memory                    80% of available RAM in GB for --metamaps [default: $params.memory]

        Workflow Options:
    --skip_host_map

        Nextflow options:
    -profile                    change the profile of nextflow both the engine and executor more details on github README
    -resume                     resume the workflow where it stopped
    """
}
