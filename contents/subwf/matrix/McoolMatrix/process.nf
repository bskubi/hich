include {withLog; stubLog} from '../../util/logs.nf'
include {buildCmd} from './functions.nf'

process MCOOL_MATRIX {
    publishDir params.general.publish.mcool ? params.general.publish.mcool : "results",
               saveAs: {params.general.publish.mcool ? it : null},
               mode: params.general.publish.mode
    label 'createMatrix'
    conda "$projectDir/env/dev_env.yml"
    container params.general.coolerContainer
    tag "$id"

    input:
    tuple val(id), val(matrixPlanName), path(pairs), path(chromsizes), val(assembly), val(matrix)

    output:
    tuple val(id), val(matrixPlanName), path(output)

    shell:
    (cmd, logMap, output) = buildCmd(id, matrixPlanName, pairs, chromsizes, assembly, matrix, task.cpus)
    withLog(cmd, logMap)

    stub:
    (cmd, logMap, output) = buildCmd(id, matrixPlanName, pairs, chromsizes, assembly, matrix, task.cpus)
    stub = "touch '${output}'"
    stubLog(stub, cmd, logMap)
}
