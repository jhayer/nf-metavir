#!/usr/bin/env nextflow

process prep_bt2_index {
    label 'prep_bt2_index'
    tag "$genome_fasta"
    publishDir "${params.output}/host_genome_idx", mode: 'copy'
    input:
        path(genome_fasta)
    output:
        path ('*.bt2')
    script:
        """
        bowtie2-build $genome_fasta ${genome_fasta.baseName}
        """
}

process bowtie2 {
    label 'bowtie2'
    publishDir "${params.output}/${id}/host_mapping", mode: 'copy'
    input:
        tuple val(id), path(illumina_clean)
        path(bt2_index_files)
    output:
        tuple val(id), path("*_un-conc.?.fastq")
        path("*_hostmap.log")
    script:
        
        index_basename = bt2_index_files[0].toString() - ~/.rev.\d.bt2?/ - ~/.\d.bt2?/
        """
        bowtie2 -x ${index_basename} -1 ${illumina_clean[0]} -2 ${illumina_clean[1]} \
            -p 4 -S ${id}_hostmap_bt2.sam --un ${id}_hostmap_un.fastq \
            --un-conc ${id}_hostmap_un-conc.fastq
        
        cp ".command.log" ${id}_hostmap.log
        """
}
