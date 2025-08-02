include {withLog; stubLog} from '../../util/logs.nf'
include {buildCmdPairtoolsParse2} from './parseHelpers.nf'

process PairtoolsParse2 {
    publishDir params.general.publish.parse ? params.general.publish.parse : "results",
               saveAs: {params.general.publish.parse ? it : null},
               mode: params.general.publish.mode

    label 'pairs'
    tag "$id"
    conda "$projectDir/env/dev_env.yml"
    container params.general.hichContainer

    input:
    tuple val(id), path(sambam), path(chromsizes), val(assembly), val(parseParams), val(sql), val(minMapq)

    output:
    tuple val(id), path("${id}.pairs.gz")

    shell:
    (cmd, logMap) = buildCmdPairtoolsParse2(id, sambam, chromsizes, assembly, parseParams, sql, minMapq, task.memory, task.cpus)
    withLog(cmd, logMap)

    stub:
    stub = "touch '${id}.pairs.gz'"
    (cmd, logMap) = buildCmdPairtoolsParse2(id, sambam, chromsizes, assembly, parseParams, sql, minMapq, task.memory, task.cpus)
    stubLog(stub, cmd, logMap)
}