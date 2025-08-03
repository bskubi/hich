include {withLog; stubLog} from '../../util/logs.nf'
include {buildCmd} from './functions.nf'

process PARSE_TO_PAIRS {
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
    tuple val(id), path(output)

    shell:
    (cmd, logMap, output) = buildCmd(id, sambam, chromsizes, assembly, parseParams, sql, minMapq, task.memory, task.cpus)
    withLog(cmd, logMap)

    stub:
    
    (cmd, logMap, output) = buildCmd(id, sambam, chromsizes, assembly, parseParams, sql, minMapq, task.memory, task.cpus)
    stub = "touch '${output}'"
    stubLog(stub, cmd, logMap)
}