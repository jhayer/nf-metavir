process megahit {
    label 'megahit'
    publishDir "${params.output}/${id}/assembly", mode: 'copy'
    input:
        tuple val(id), path(illumina)
    output:
        tuple val(id), path("megahit/final.contigs.fa")
        path("megahit")
    script:
        """
        megahit -t ${task.cpus} -m 4000 -1 ${illumina[0]} \
        -2 ${illumina[1]} -o megahit
        """
}
