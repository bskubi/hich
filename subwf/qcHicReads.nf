process PairtoolsStats {
    publishDir params.general.publish.pair_stats ? params.general.publish.pair_stats : "results",
               saveAs: {params.general.publish.pair_stats ? it : null},
               mode: params.general.publish.mode
    container "bskubi/hich:latest"

    input:
    tuple val(id), val(pairs_id), path(pairs)

    output:
    tuple val(id), val(pairs_id), path("${id}.${pairs_id}.stats.txt")

    shell:
    "pairtools stats --output ${id}.${pairs_id}.stats.txt ${pairs}"

    stub:
    "touch ${id}.${pairs_id}.stats.txt"
}

process MultiQC {
    publishDir params.general.publish.qc ? params.general.publish.qc : "results",
               saveAs: {params.general.publish.qc ? it : null},
               mode: params.general.publish.mode
    container "bskubi/hich:latest"

    input:
    tuple val(report_name), path(stats)

    output:
    path("${report_name}.multiqc_report.html")

    shell:
    "multiqc --force --filename ${report_name}.multiqc_report.html --module pairtools . && chmod -R a+w ."

    stub:
    "touch ${report_name}.multiqc_report.html"
}

workflow QCReads {
    take:
    samples
    report_name
    
    main:
    samples
    | map{
        sample ->
        pairs = ["pairs",
                  "frag_pairs",
                  "dedup_pairs",
                  "select_pairs",
                  "biorep_merge_pairs",
                  "condition_merge_pairs"]
        
        emit = pairs.collect{
            pairfile ->
            
            sample.containsKey(pairfile) && sample.get(pairfile).getClass() == sun.nio.fs.UnixPath ? [sample.id, pairfile, sample[pairfile]] : null
        }.findAll{it != null}
        emit
    }
    | collect
    | flatMap
    | PairtoolsStats
    | map{it[2]}
    | collect
    | map{[report_name, it]}
    | MultiQC
    

    emit:
    samples

}