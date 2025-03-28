include {withLog; stubLog; skip} from '../extraops.nf'


process PairtoolsStats {
    publishDir params.general.publish.pairStats ? params.general.publish.pairStats : "results",
               saveAs: {params.general.publish.pairStats ? it : null},
               mode: params.general.publish.mode
    conda "bioconda::pairtools"
    container params.general.hichContainer
    label 'pairs'
    cpus 8

    input:
    tuple val(id), val(pairs_id), path(pairs)

    output:
    tuple val(id), val(pairs_id), path("${id}.${pairs_id}.stats.txt")

    shell:
    cmd = "pairtools stats --output '${id}.${pairs_id}.stats.txt'  --nproc-in ${task.cpus} --nproc-out ${task.cpus} '${pairs}'"

    logMap = [
        task: "PairtoolsStats",
        input: [
            id: id,
            pairs_id: pairs_id,
            pairs: pairs
        ],
        output: [
            stats: "${id}.${pairs_id}.stats.txt"
        ]
    ]

    withLog(cmd, logMap)

    stub:
    stub = "touch '${id}.${pairs_id}.stats.txt'"

    cmd = "pairtools stats --output '${id}.${pairs_id}.stats.txt'  --nproc-in ${task.cpus} --nproc-out ${task.cpus} '${pairs}'"

    logMap = [
        task: "PairtoolsStats",
        input: [
            id: id,
            pairs_id: pairs_id,
            pairs: pairs
        ],
        output: [
            stats: "${id}.${pairs_id}.stats.txt"
        ]
    ]

    stubLog(stub, cmd, logMap)
}

process MultiQC {
    publishDir params.general.publish.qc ? params.general.publish.qc : "results",
               saveAs: {params.general.publish.qc ? it : null},
               mode: params.general.publish.mode
    conda "env/multiqc.yml"
    container "bskubi/hich:latest"

    input:
    tuple val(report_name), path(stats)

    output:
    path("${report_name}.multiqc_report.html")

    shell:
    cmd = "multiqc --force --no-version-check --filename '${report_name}.multiqc_report.html' --module pairtools . && chmod -R a+w ."
    cmd

    stub:
    stub = "touch '${report_name}.multiqc_report.html'"
    cmd = "multiqc --force --no-version-check --filename '${report_name}.multiqc_report.html' --module pairtools . && chmod -R a+w ."
    logMap = [
        task: task,
        input: [
            report_name: report_name,
            stats: stats
        ],
        output: [
            report: "${report_name}.multiqc_report.html"
        ]
    ]

    stubLog(stub, cmd, logMap)
}

workflow QCReads {
    take:
    samples
    report_name
    
    main:
    samples
    | filter{!skip("qcReads")}
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