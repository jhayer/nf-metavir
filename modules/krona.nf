process krona_chart_kraken {
    label 'krona_chart_kraken'
    publishDir "${params.output}/${id}/taxonomic_classif", mode: 'copy'
    input:
        path(krak_report)
    output:
        path("*krona.html")
    script:
        """
        suffix="report.txt"
        prefix=${krak_report%$suffix}
        ## parse kraken-report
        parse_to_krona_v2.py ${krak_report}
        ## create kron files
        cut -d$'\t' -f3,6- "${prefix}report.parsed.txt" > "${prefix}krona.in"
        ktImportText -o "${prefix}krona.html" "${prefix}krona.in"

        rm "${prefix}report.parsed.txt"
        rm "${prefix}krona.in"
        """
}
