include {withLog; stubLog} from '../../util/logs.nf'
include {buildCmd} from './functions.nf'

process SELECT_PAIRS {
    publishDir params.general.publish.select ? params.general.publish.select : "results",
               saveAs: {params.general.publish.select ? it : null},
               mode: params.general.publish.mode

    label 'pairs'
    tag "$id"
    conda "$projectDir/env/dev_env.yml"
    container params.general.hichContainer

    input:
    tuple val(id), path(pairs), val(selectPairs_opts)

    output:
    tuple val(id), path(output)

    shell:
    (cmd, logMap, output) = buildCmd(id, pairs, selectPairs_opts, task.cpus)
    withLog(cmd, logMap)

    stub:
    
    (cmd, logMap, output) = buildCmd(id, pairs, selectPairs_opts, task.cpus)
    stub = "touch '${output}'"
    stubLog(stub, cmd, logMap)
}