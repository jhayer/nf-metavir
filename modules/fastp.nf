process fastp {
    label 'fastp'
    publishDir "${params.output}/${name}/qc", mode: 'copy', pattern: "*_R*_clean.fastq"
    input:
        tuple val(name), path(illumina)
    output:
        tuple val(name), path("*_R?_clean.fastq")
    script:
        """
        fastp -i ${illumina[0]} -I ${illumina[1]} -o ${name}_R1_clean.fastq -O ${name}_R2_clean.fastq
        """
}
