include {withLog; stubLog} from '../../util/logs.nf'
include {buildCmd} from './functions.nf'

process INGEST_PAIRS {
    publishDir params.general.publish.flip_sort ? params.general.publish.flip_sort : "results",
               saveAs: {params.general.publish.flip_sort ? it : null},
               mode: params.general.publish.mode

    conda "$projectDir/env/dev_env.yml"
    container params.general.hichContainer
    label 'pairs'
    tag "$id"


    input:
    tuple val(id), path(pairs), path(chromsizes), val(ingest_pairs_opts)

    output:
    tuple val(id), path(output)

    shell:
    (cmd, logMap, output) = buildCmd(id, pairs, chromsizes, ingest_pairs_opts, task.cpus, task.memory)
    withLog(cmd, logMap)

    stub:
    (cmd, logMap, output) = buildCmd(id, pairs, chromsizes, ingest_pairs_opts, task.cpus, task.memory)
    stub = "touch '${output}'"
    stubLog(stub, cmd, logMap)
}
