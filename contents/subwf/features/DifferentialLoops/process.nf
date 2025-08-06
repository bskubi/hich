include {withLog; stubLog} from '../../util/logs.nf'
include {buildCmd} from './functions.nf'

process DIFFERENTIAL_LOOPS {
    container params.general.mustacheContainer
    label 'features'
    tag "$id1 $id2"

    input:
    tuple val(id1), val(id2), val(prefix), path(matrix1, stageAs: '__matrix1__/*'), path(matrix2, stageAs: '__matrix1__/*'), val(analysisPlan)

    output:
    tuple path(loop1), path(loop2), path(diffloop1), path(diffloop2)

    shell:
    (cmd, logMap, loop1, loop2, diffloop1, diffloop2) = buildCmd(id1, id2, prefix, matrix1, matrix2, analysisPlan)
    withLog(cmd, logMap)

    stub:
    (cmd, logMap, loop1, loop2, diffloop1, diffloop2) = buildCmd(id1, id2, prefix, matrix1, matrix2, analysisPlan)
    stub = "touch '${loop1}' '${loop2}' '${diffloop1}' '${diffloop2}'"
    stubLog(stub, cmd, logMap)
}