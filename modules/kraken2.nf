process kraken2prot_reads {
    label 'kraken2prot_reads'
    publishDir "${params.output}/${id}/taxonomic_classif/reads", mode: 'copy'
    input:
        tuple val(id), path(illumina)
    output:
        path("kn2_prot_nr")
    script:
        """
        kraken2 --db ${params.k2prot_db} --memory-mapping \
            --threads ${task.cpus} --output ${id}_kn2_nr-res.txt \
            --report ${id}_kn2_nr-report.txt \
            --paired ${illumina[0]} ${illumina[1]}
        """
}
