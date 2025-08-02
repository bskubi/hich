include {withLog; stubLog} from '../../util/logs.nf'
include {buildCmdPairtoolsSelect} from './selectHelpers.nf'

process PairtoolsSelect {
    publishDir params.general.publish.select ? params.general.publish.select : "results",
               saveAs: {params.general.publish.select ? it : null},
               mode: params.general.publish.mode

    label 'pairs'
    tag "$id"
    conda "$projectDir/env/dev_env.yml"
    container params.general.hichContainer

    input:
    tuple val(id), path(pairs), val(pairtoolsSelectParams), val(pairtoolsSelectFilters)

    output:
    tuple val(id), path("${id}_select.pairs.gz")

    shell:
    (cmd, logMap) = buildCmdPairtoolsSelect(id, pairs, pairtoolsSelectParams, pairtoolsSelectFilters, task.cpus)
    withLog(cmd, logMap)

    stub:
    stub = "touch '${id}_select.pairs.gz'"
    (cmd, logMap) = buildCmdPairtoolsSelect(id, pairs, pairtoolsSelectParams, pairtoolsSelectFilters, task.cpus)
    stubLog(stub, cmd, logMap)
}