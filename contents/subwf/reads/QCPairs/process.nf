include {withLog; stubLog} from '../../util/logs.nf'

process PAIRTOOLS_STATS {
    publishDir params.general.publish.pairStats ? params.general.publish.pairStats : "results",
               saveAs: {params.general.publish.pairStats ? it : null},
               mode: params.general.publish.mode
    label 'stats'
    tag "$id"
    conda "$projectDir/env/dev_env.yml"
    container params.general.hichContainer


    input:
    tuple val(id), val(pairsID), val (statsKey), path(pairs)

    output:
    tuple val(id), val(pairsID), val (statsKey), path("${id}.${pairsID}.stats.txt")

    shell:
    def cmd = "pairtools stats --output '${id}.${pairsID}.stats.txt'  --nproc-in ${task.cpus} --nproc-out ${task.cpus} '${pairs}'"

    def logMap = [
        task: "PairtoolsStats",
        input: [
            id: id,
            pairsID: pairsID,
            pairs: pairs
        ],
        output: [
            stats: "${id}.${pairsID}.stats.txt"
        ]
    ]

    withLog(cmd, logMap)

    stub:
    def stub = "touch '${id}.${pairsID}.stats.txt'"

    def cmd = "pairtools stats --output '${id}.${pairsID}.stats.txt'  --nproc-in ${task.cpus} --nproc-out ${task.cpus} '${pairs}'"

    def logMap = [
        task: "PairtoolsStats",
        input: [
            id: id,
            pairsID: pairsID,
            pairs: pairs
        ],
        output: [
            stats: "${id}.${pairsID}.stats.txt"
        ]
    ]

    stubLog(stub, cmd, logMap)
}

process MULTIQC_PAIRTOOLS {
    publishDir params.general.publish.qc ? params.general.publish.qc : "results",
               saveAs: {params.general.publish.qc ? it : null},
               mode: params.general.publish.mode
    label 'stats'
    tag "$reportName"
    conda "$projectDir/env/dev_env.yml"
    container params.general.multiqcContainer
    

    input:
    tuple val(reportName), path(stats)

    output:
    path("${reportName}.multiqcReport.html")

    shell:
    def cmd = "multiqc --force --no-version-check --filename '${reportName}.multiqcReport.html' --module pairtools . && chmod -R a+w ."
    def logMap = [
        task: task,
        input: [
            reportName: reportName,
            stats: stats
        ],
        output: [
            report: "${reportName}.multiqcReport.html"
        ]
    ]

    withLog(cmd, logMap)

    stub:
    def stub = "touch '${reportName}.multiqcReport.html'"
    def cmd = "multiqc --force --no-version-check --filename '${reportName}.multiqcReport.html' --module pairtools . && chmod -R a+w ."
    def logMap = [
        task: task,
        input: [
            reportName: reportName,
            stats: stats
        ],
        output: [
            report: "${reportName}.multiqcReport.html"
        ]
    ]

    stubLog(stub, cmd, logMap)
}