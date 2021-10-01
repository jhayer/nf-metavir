process diamond_contigs{
    label 'diamond_contigs'
    publishDir "${params.output}/${id}/taxonomic_classif/contigs", mode: 'copy'
    input:
        tuple val(id), path(contigs)
        path(db_diamond)
        path(kraken1_nt_db)
        path(kraken_report)
    output:
        path("${id}_dx_tax.tab")
        path("${id}_dx_krak-report.txt")
    script:
        """
        #Use of format 102 for producing a taxonomy table output file (it uses LCA)
        #This output will later be made compatible for kraken-report for import into Pavian
        diamond blastx -d ${db_diamond} -q ${contigs} -o ${id}_dx_tax.tab -f 102 -p 8

        # adding the first column U/C to the tab file for kraken-report
        awk -F'\\t' '{if(\$2>0)\$1="C" FS \$1;else \$1="U" FS \$1;}1' OFS='\\t' ${id}_dx_tax.tab > ${id}_dx_tax_UC.tab

        #the diamond output can now be converted into kraken report for pavian
        ./${kraken_report} --db ${kraken1_nt_db} ${id}_dx_tax_UC.tab > ${id}_dx_krak-report.txt
        rm ${id}_dx_tax_UC.tab
        """
}

process diamond4megan_contigs{
    label 'diamond4megan_contigs'
    publishDir "${params.output}/${id}/taxonomic_classif/contigs", mode: 'copy'
    input:
        tuple val(id), path(contigs)
        path(db_diamond)
    output:
        path("${id}_dx.daa")
    script:
        """
        #Use of format 100 for producing a daa file compatible with Megan
        diamond blastx -d ${db_diamond} -q ${contigs} -o ${id}_dx.daa -f 100 -p ${task.cpus}
        """
}
