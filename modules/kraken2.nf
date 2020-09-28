process kraken2prot_reads {
    label 'kraken2prot_reads'
    publishDir "${params.output}/${id}/taxonomic_classif/reads", mode: 'copy'
    input:
        tuple val(id), path(illumina)
        path(db_k2prot)
    output:
        path("kn2_prot_nr")
    script:
        """
        kraken2 --db ${db_k2prot} --memory-mapping \
            --threads ${task.cpus} --output ${id}_kn2_nr-res.txt \
            --report ${id}_kn2_nr-report.txt \
            --paired ${illumina[0]} ${illumina[1]}
        """
}
