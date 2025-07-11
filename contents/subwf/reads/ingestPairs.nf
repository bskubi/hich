include {QCPairs} from './qcPairs.nf'
include {emptyOnLastStep; skip} from '../util/cli.nf'
include {keyUpdate} from '../util/keyUpdate.nf'
include {withLog; stubLog} from '../util/logs.nf'

process PairtoolsFlipSort {
    publishDir params.general.publish.flip_sort ? params.general.publish.flip_sort : "results",
               saveAs: {params.general.publish.flip_sort ? it : null},
               mode: params.general.publish.mode

    label 'pairs'
    conda "$projectDir/env/dev_env.yml"

    input:
    tuple val(id), path(pairs), path(chromsizes), val(sql)

    output:
    tuple val(id), path("${id}.pairs.gz")

    shell:

    reshapeCmd = sql ? ["hich pairs sql --memory-limit '${task.memory}' --threads '${task.cpus}' '${sql}' '${pairs}'"] : []
    flipCmd = ["pairtools flip --chroms-path '${chromsizes}'  --nproc-in ${task.cpus} --nproc-out ${task.cpus}"]
    
    memory = task.memory ? task.memory.toGiga() : "2"
    sortCmd = ["pairtools sort --output '${id}.pairs.gz'  --memory ${memory}G --nproc-in ${task.cpus} --nproc-out ${task.cpus}"]

    if (!reshapeCmd) {
        flipCmd = flipCmd + ["'${pairs}'"]
        flipCmd = [flipCmd.join(" ")]
    }

    cmdParts = reshapeCmd + flipCmd + sortCmd
    cmd = cmdParts.join(" | ")

    logMap = [task: "PairtoolsFlipSort", input: [id: id, pairs: pairs, chromsizes: chromsizes, sql: sql], 
    output: [pairs: "${id}.pairs.gz"]]
    withLog(cmd, logMap)

    stub:
    stub = "touch ${id}.pairs.gz"

    reshapeCmd = sql ? ["hich pairs sql --memory-limit '${task.memory}' --threads '${task.cpus}' '${sql}' '${pairs}'"] : []
    flipCmd = ["pairtools flip --chroms-path '${chromsizes}'  --nproc-in ${task.cpus} --nproc-out ${task.cpus}"]
    memory = task.memory ? task.memory.toGiga() : "2"
    sortCmd = ["pairtools sort --output '${id}.pairs.gz'  --memory ${memory}G --nproc-in ${task.cpus} --nproc-out ${task.cpus}"]

    if (!reshapeCmd) {
        flipCmd = flipCmd + ["'${pairs}'"]
        flipCmd = [flipCmd.join(" ")]
    }

    cmdParts = reshapeCmd + flipCmd + sortCmd
    cmd = cmdParts.join(" | ")

    logMap = [task: "PairtoolsFlipSort", input: [id: id, pairs: pairs, chromsizes: chromsizes, sql: sql], 
    output: [pairs: "${id}.pairs.gz"]]
    stubLog(stub, cmd, logMap)
}

workflow IngestPairs {
    take:
        samples

    main:


    samples
        | filter{!skip("ingestPairs") && it.datatype == "pairs"}
        | map{tuple(it.id, it.pairs, it.chromsizes, it.reshapeParams)}
        | PairtoolsFlipSort
        | map{[id:it[0], pairs:it[1], latest:it[1], latestPairs:it[1]]}
        | set{result}
    keyUpdate(samples, result, "id") | set{samples}

    if ("ingestPairs" in params.general.get("qcAfter")) {
        samples = QCPairs(samples, ["pairs"], "ingestPairs")
    }

    samples = emptyOnLastStep("ingestPairs", samples)

    emit:
        samples
}