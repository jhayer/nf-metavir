process kraken2prot_reads {
    label 'kraken2prot_reads'
    publishDir "${params.output}/${id}/taxonomic_classif/reads", mode: 'copy'
    input:
        tuple val(id), path(illumina)
        path(db_k2prot)
    output:
        path("*_kn2_nr-re*.txt")
    script:
        """
        kraken2 --db ${db_k2prot} --memory-mapping \
            --threads ${task.cpus} --output ${id}_kn2_nr-res.txt \
            --report ${id}_kn2_nr-report.txt \
            --paired ${illumina[0]} ${illumina[1]}
        """
}

process kraken2prot_contigs {
    label 'kraken2prot_contigs'
    publishDir "${params.output}/${id}/taxonomic_classif/contigs", mode: 'copy'
    input:
        tuple val(id), path(contigs)
        path(db_k2prot)
    output:
        path("*_kn2_nr-re*.txt")
    script:
        """
        kraken2 --db ${db_k2prot} --memory-mapping \
            --threads ${task.cpus} --output ${id}_kn2_nr-res.txt \
            --report ${id}_kn2_nr-report.txt ${contigs}
        """
}
