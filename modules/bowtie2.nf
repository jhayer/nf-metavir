process prep_bt2_index {
    label 'prep_bt2_index'
    publishDir "${params.output}/host_genome_idx/", mode: 'copy'
    input:
        val(fasta)
    output:
        set val("bowtieIndex"), file("*") into bowtie_index
    script:
        """
        bowtie2-build ${fasta} bowtie_index
        """
}

process bowtie2 {
    label 'bowtie2'
    publishDir "${params.output}/${id}/host_mapping", mode: 'copy'
    input:
        tuple val(id), path(illumina_clean)
    output:
        tuple val(id), path("*_R?_un-conc.fastq")
        path("${id}_bt2.log")
    script:
        """
        bowtie2 -x ${bowtie_index} -1 ${illumina_clean[0]} -2 ${illumina_clean[1]} \
            -p 4 -S ${id}_hostmap_bt2.sam --un ${id}_hostmap_un.fastq \
            --un-conc ${id}_hostmap_un-conc.fastq --met-file ${id}_bt2.log
        """
}
