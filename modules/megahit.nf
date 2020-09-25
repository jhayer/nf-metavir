process megahit {
    label 'megahit'
    publishDir "${params.output}/${id}", mode: 'copy'
    input:
        tuple val(id), path(illumina)
    output:
        tuple val(id), path("${id}_megahit/final.contigs.fa")
        path("${id}_megahit")
    script:
        """
        megahit -t ${task.cpus} -m ${task.memory} -1 ${illumina[0]} \
        -2 ${illumina[1]} -o "${id}_megahit"
        """
}
