include {skip} from '../util/cli.nf'
include {withLog; stubLog} from '../util/logs.nf'
include {keyUpdate} from '../util/keyUpdate.nf'

process PairtoolsStats {
    publishDir params.general.publish.pairStats ? params.general.publish.pairStats : "results",
               saveAs: {params.general.publish.pairStats ? it : null},
               mode: params.general.publish.mode
    label 'stats'
    tag "$id"

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

process MultiQC {
    publishDir params.general.publish.qc ? params.general.publish.qc : "results",
               saveAs: {params.general.publish.qc ? it : null},
               mode: params.general.publish.mode
    label 'stats'
    tag "$reportName"

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

workflow QCPairs {
    take:
    samples
    pairs
    reportName
    
    main:

    if (!skip("qcpairs")) {
        samples
            | map{
                sample ->
                
                pairs.collect{
                    pairfile ->
                    hasPairFile = sample.containsKey(pairfile)
                    pairFileIsPath = sample.get(pairfile).getClass() == sun.nio.fs.UnixPath
                    statsKey = pairfile + "Stats"
                    noStatsFile = !sample.containsKey(statsKey)
                    computeStats = hasPairFile && pairFileIsPath && noStatsFile
                    statsInput = [sample.id, pairfile, statsKey, sample[pairfile]]
                    computeStats ? statsInput : null
                }.findAll{it != null}
            }
            | collect
            | flatMap
            | PairtoolsStats
            | map{it[3]}
            | collect
            | map{paths -> [reportName, paths]}
            | MultiQC
    }



    emit:
    samples

}