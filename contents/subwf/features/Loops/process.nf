include {withLog; stubLog} from '../../util/logs.nf'
include {buildCmd} from './functions.nf'

process LOOPS{
    publishDir "results/loops",
               mode: params.general.publish.mode
    container params.general.mustacheContainer
    label 'features'
    tag "$id"

    input:
    tuple val(id), path(matrix), val(analysisPlan)

    output:
    path(output)

    shell:
    (cmd, logMap, output) = buildCmd(id, matrix, analysisPlan)
    withLog(cmd, logMap)

    stub:
    (cmd, logMap, output) = buildCmd(id, matrix, analysisPlan)
    stub = "touch '${output}'"
    stubLog(stub, cmd, logMap)
}