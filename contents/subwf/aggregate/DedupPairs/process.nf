include {withLog; stubLog} from '../../util/logs.nf'
include {buildCmd} from './functions.nf'

process DEDUP_PAIRS {
    publishDir params.general.publish.parse ? params.general.publish.dedup : "results",
               saveAs: {params.general.publish.dedup ? it : null},
               mode: params.general.publish.mode

    label 'pairs'
    conda "$projectDir/env/dev_env.yml"
    container params.general.hichContainer

    input:
    tuple val(id), path(pairs), val(singleCell), val(maxMismatch), val(method), val(pairtoolsDedupParams)

    output:
    tuple val(id), path(output)

    shell:
    (cmd, logMap, output) = buildCmd(id, pairs, singleCell, maxMismatch, method, pairtoolsDedupParams, task.cpus)
    withLog(cmd, logMap)

    stub:
    (cmd, logMap, output) = buildCmd(id, pairs, singleCell, maxMismatch, method, pairtoolsDedupParams, task.cpus)
    stub = "touch '${output}'"
    stubLog(stub, cmd, logMap)
}
