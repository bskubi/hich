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
    tuple val(id), val(matrixPlanName), path(pairs), path(chromsizes), val(assembly), val(matrix), val(coolerCloadParams), val(coolerZoomifyParams)

    output:
    tuple val(id), val(matrixPlanName), path("${id}.mcool")

    shell:
    (cmd, logMap, output) = buildCmd(id, matrixPlanName, pairs, chromsizes, assembly, matrix, coolerCloadParams, coolerZoomifyParams, task.cpus)
    withLog(cmd, logMap)

    stub:
    (cmd, logMap, output) = buildCmd(id, matrixPlanName, pairs, chromsizes, assembly, matrix, coolerCloadParams, coolerZoomifyParams, task.cpus)
    stub = "touch '${output}'"
    stubLog(stub, cmd, logMap)
}
