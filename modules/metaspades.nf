process metaspades {
    label 'metaspades'
    publishDir "${params.output}/${id}/assembly", mode: 'copy'
    input:
        tuple val(id), path(illumina)
    output:
        tuple val(id), path("metaspades/scaffolds.fasta")
        path("metaspades")
    script:
        """
        metaspades.py -1 ${illumina[0]} -2 ${illumina[1]} \
         -o metaspades -t ${task.cpus}
        """
}
