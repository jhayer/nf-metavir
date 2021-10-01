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
    nextflow run main.nf --illumina illumina/ -profile planet --skip_host_map --k2prot_db k2nr_path
    or
    nextflow run main.nf --illumina illumina/ --host_map host_ref/ -profile uppmax
        Input:
    --illumina                  path to the directory containing the illumina read file (fastq) (default: $params.illumina)
        Optional input:
    --host_ref                  path to the host reference genome to map on [default: $params.host_ref]
    --k2nt_db                   path to the Kraken2 nucleotide database (e.g. nt) [default: $params.k2nt_db]
    --k2prot_db                 path to the Kraken2 protein database (e.g. nr) [default: $params.k2prot_db]
    --diamond_db                path to the Diamond protein database (e.g. nr)
    --kraken1report             path to Kraken 1 program kraken-report [default: kraken-report]

        Output:
    --output                    path to the output directory (default: $params.output)

        Outputed directories:
    qc                          The reads file after qc, qc logs and host mapping logs
    assembly                    The megahit assembly output directory
    taxonomic_classif
        reads                       The taxonomic classifications at reads level
        contigs                     The taxonomic classifications at contigs level

        Basic Parameter:
    --cpus                      max cores for local use [default: $params.cpus]
    --memory                    80% of available RAM in GB [default: $params.memory]

        Workflow Options:
    --skip_host_map             if set, no host mapping.[default: $params.skip_host_map]
    --diamond4megan             produce a diamond daa file for Megan (need to set diamond_db)[default: $params.diamond4megan]
    --skip_diamond4pavian       skip the run of Diamond for Pavian output [default: $params.skip_diamond4pavian]

        Nextflow options:
    -profile                    change the profile of nextflow both the engine and executor more details on github README
    -resume                     resume the workflow where it stopped
    """
}

workflow {

    // error handling
    if (
        workflow.profile.contains('planet') ||
        workflow.profile.contains('uppmax') ||
        workflow.profile.contains('itrop') ||
        workflow.profile.contains('local')
    ) { "executer selected" }
    else { exit 1, "No executer selected: -profile local/uppmax/planet/itrop"}


    //*************************************************
    // STEP 0 - Include needed modules
    //*************************************************

    include {fastp} from './modules/fastp' params(output: params.output)
    // including bowtie2 module if host mapping needed
    if (params.skip_host_map==false) {
        // do the host mapping
        include {prep_bt2_index} from './modules/bowtie2' params(output: params.output)
        include {bowtie2} from './modules/bowtie2' params(output: params.output)
    }
    else {
            if(params.host_ref) {
                exit 1, "skip_host_map options and host_ref are incompatible"
            }
    }
    // including megahit module
    include {megahit} from './modules/megahit' params(output: params.output)

    // including Kraken2 - protein level
    if (params.k2prot_db){
        include {kraken2prot_reads} from './modules/kraken2.nf' params(output: params.output)
        include {kraken2prot_contigs} from './modules/kraken2.nf' params(output: params.output)
    }
    // including Kraken2 - nucleotide level
    if (params.k2nt_db) {
        include {kraken2nt_reads} from './modules/kraken2.nf' params(output: params.output)
        include {kraken2nt_contigs} from './modules/kraken2.nf' params(output: params.output)
    }
    // TODO: here add warnings if no K2 db selected

    if (params.krona_chart_kraken) {
        include {krona_chart_kraken} from './modules/krona.nf' params(output: params.output)
    }

    // including Diamond: ouputs for Pavian or Megan if needed
    if (params.diamond_db) {
        if (params.skip_diamond4pavian==false){
            include {diamond_contigs} from './modules/diamond.nf' params(output: params.output)
        }

        if (params.diamond4megan==true) {
            include {diamond4megan_contigs} from './modules/diamond.nf' params(output: params.output)
        }
    }
    else {
        if(params.diamond4megan==true){
            exit 1, "You need to specify a Diamond database to use"
        }
    }

    //*************************************************
    // STEP 1 QC with fastp
    //*************************************************

    // DATA INPUT ILLUMINA
    illumina_input_ch = Channel
        .fromFilePairs( "${params.illumina}/*_{1,2}.fastq{,.gz}", checkIfExists: true)
        .view()

    // run fastp module
    fastp(illumina_input_ch)
    illumina_clean_ch = fastp.out[0]

    //*************************************************
    // STEP 2 - Optional - Host mapping
    //*************************************************
    if (params.skip_host_map==false){
            if(params.host_ref) {
                fasta = file(params.host_ref)
                if( !fasta.exists() ) exit 1, "Fasta file not found: ${params.host_ref}"

                // The reference genome file
                Channel.fromPath(params.host_ref, checkIfExists: true)
                        .ifEmpty { exit 1, "Cannot find host genome matching ${params.host_ref}!\n" }
                        .set {genome}

                if (
                file("${fasta.baseName}.1.bt2").exists() &&
                file("${fasta.baseName}.2.bt2").exists() &&
                file("${fasta.baseName}.3.bt2").exists() &&
                file("${fasta.baseName}.4.bt2").exists() &&
                file("${fasta.baseName}.rev.1.bt2").exists() &&
                file("${fasta.baseName}.rev.2.bt2").exists()
                ){
                    index_files = Channel.fromPath("${fasta.baseName}*.bt2")
                }
                else {
                    //prep genome index
                    prep_bt2_index(genome)
                    index_files = prep_bt2_index.out
                    //bowtie2(illumina_clean_ch,prep_bt2_index.out.collect())
                }
                // mapping bowtie2
                bowtie2(illumina_clean_ch,index_files.collect())
                illumina_host_unmapped_ch = bowtie2.out[0]
            }
            else {
                "Warning: No reference genome specified! Skipping host mapping"
                illumina_host_unmapped_ch = illumina_clean_ch
            }
    }
    else {
        illumina_host_unmapped_ch = illumina_clean_ch
    }

    //*************************************************
    // STEP 3 - Assembly
    //*************************************************
    megahit(illumina_host_unmapped_ch)
    contigs_ch = megahit.out[0]

    //*************************************************
    // STEP 4A - taxonomic classification prot level
    //           on reads and contigs - Kraken2
    //*************************************************
    if (params.k2prot_db){
        db_k2prot = file(params.k2prot_db)
        kraken2prot_reads(illumina_host_unmapped_ch, db_k2prot)
        kraken2prot_contigs(contigs_ch, db_k2prot)
        if (params.krona_chart==true) {
            k2res_reads_rep_ch = kraken2prot_reads.out[1]
        //    krona_chart_kraken(k2res_reads_rep_ch)
        }
    }

    //*************************************************
    // STEP 4B - taxonomic classification nucleotide
    //           on reads and contigs - Kraken2
    //*************************************************
    if (params.k2nt_db) {
        db_k2nt = file(params.k2nt_db)
        kraken2nt_reads(illumina_host_unmapped_ch, db_k2nt)
        kraken2nt_contigs(contigs_ch, db_k2nt)
    }

    //*************************************************
    // STEP 4C - taxonomic classification protein level
    //           on reads and contigs - diamond
    //*************************************************
    if (params.diamond_db) {
        db_diamond = file(params.diamond_db)
        // on contigs only for now as Diamond does not support paired ends
        if (params.skip_diamond4pavian==false){
            // run diamond with output compatible for pavian
            //need kraken and kraken db for kraken_report
            kraken1_nt_db = file(params.krak1_nt_db)
            kraken_report = file(params.kraken1report)
            diamond_contigs(contigs_ch, db_diamond, kraken1_nt_db, kraken_report)
        }
        // run diamond with daa output compatible for Megan
        if (params.diamond4megan==true) {
            diamond4megan_contigs(contigs_ch, db_diamond)
        }
    }

    //*************************************************
    // STEP 4D - taxonomic classification protein level
    //           on reads and contigs - Kaiju
    //*************************************************


    //*************************************************
    // STEP 5 - taxonomic visualisation - KronaChart
    //*************************************************


}
