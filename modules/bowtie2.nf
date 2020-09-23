process bowtie2 {
    label 'bowtie2'
    publishDir "${params.output}/${id}/host_mapping", mode: 'copy'
    input:
        tuple val(id), path(illumina_clean), path(ref_genome)
        file fasta from ref_genome
    output:
        tuple val(id), path("*_R?_un-conc.fastq")
        path("${id}_bt2.log")
    script:
        """
        bowtie2-build ${fasta} ref_index

        bowtie2 -x ${ref_index} -1 ${illumina_clean[0]} -2 ${illumina_clean[1]} \
            -p 4 -S ${id}_hostmap_bt2.sam --un ${id}_hostmap_un.fastq \
            --un-conc ${id}_hostmap_un-conc.fastq --met-file ${id}_bt2.log
        """
}
