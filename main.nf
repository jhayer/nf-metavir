#!/usr/bin/env nextflow
nextflow.enable.dsl=2

start_var = Channel.from("""
********* Start running Metavir pipeline *********
nf-metavir is a workflow for metagenomics qc, assembly and taxonomic classification
**************************************
""")
start_var.view()

if (params.help) { exit 0, helpMSG() }

// Help Message
def helpMSG() {
    log.info """
    ********* Assembly and taxonomic classification workflow for (viral) metagenomics *********

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

if( !nextflow.version.matches('20.+') ) {
    ch_ver=Channel.from("This workflow requires Nextflow version 20.07 or greater -- You are running version $nextflow.version").view()
    exit 1
}

workflow {

    // error handling
    if (
        workflow.profile.contains('planet') ||
        workflow.profile.contains('uppmax')
    ) { "executer selected" }
    else { exit 1, "No executer selected: -profile uppmax or -profile planet"}


    include {fastp} from './modules/fastp' params(output: params.output)

    if (params.skip_host_map==false) {
        // do the host mapping
        include {prep_bt2_index} from './modules/bowtie2' params(output: params.output)
        include {bowtie2} from './modules/bowtie2' params(output: params.output)
    }


    //*************************************************
    // STEP 1 QC with fastp
    //*************************************************

    // DATA INPUT ILLUMINA
    illumina_input_ch = Channel
        .fromFilePairs( "${params.illumina}/*_R{1,2}_001.fastq{,.gz}", checkIfExists: true)
        .view()

    // run fastp module
//    fastp(illumina_input_ch)
//    illumina_clean_ch = fastp.out

    //*************************************************
    // STEP 2 - Optional - Host mapping
    //*************************************************
    if (params.skip_host_map==false){
            if(params.host_ref) {
                fasta = file(params.host_ref)
                if( !fasta.exists() ) exit 1, "Fasta file not found: ${params.host_ref}"
            }
            else {
                exit 1, "No reference genome specified! A fasta file is required"
            }

            // The reference genome file
           // genome_file = file(params.host_ref)
           // genome_file = Channel.fromPath(params.host_ref)
            Channel.fromPath(params.host_ref, checkIfExists: true)
                    .ifEmpty { exit 1, "Cannot find host genome matching ${params.host_ref}!\n" }
                    .set {genome}
            //genome_file.view()
            bt2_index(genome)
          //  bowtie2(illumina_clean_ch)
          //  illumina_host_unmapped_ch = bowtie2.out
    }

}

process bt2_index {
    tag "$genome_fasta"
    publishDir "${params.output}/host_genome_idx", mode: 'copy'
    input:
        path(genome_fasta)
    output:
        path ('*.bt2')
    script:
        """
        bowtie2-build $genome_fasta ${genome_fasta.baseName}
        """
}


