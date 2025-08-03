include {withLog; stubLog} from '../../util/logs.nf'
include {buildCmd} from './functions.nf'

process JuicerToolsPre {
    /*
        Juicer Tools Pre documentation: https://github.com/aidenlab/juicer/wiki/Pre
    */
    publishDir params.general.publish.hic ? params.general.publish.hic : "results",
               saveAs: {params.general.publish.hic ? it : null},
               mode: params.general.publish.mode
    tag "$id"

    label 'createMatrix'
    conda "$projectDir/env/dev_env.yml"
    container params.general.juicerContainer

    input:
    tuple val(id), val(matrixPlanName), path(pairs), path(chromsizes), val(pairsFormat), val(matrix), val(juicerToolsPreParams), val(flags)

    output:
    tuple val(id), val(matrixPlanName), path(output)

    shell:
    (cmd, logMap, output) = buildCmd(id, matrixPlanName, pairs, chromsizes, pairsFormat, matrix, juicerToolsPreParams, flags, task.memory, task.cpus)
    withLog(cmd, logMap)

    stub:
    (cmd, logMap, output) = buildCmd(id, matrixPlanName, pairs, chromsizes, pairsFormat, matrix, juicerToolsPreParams, flags, task.memory, task.cpus)
    stub = "touch '${output}'"
    stubLog(stub, cmd, logMap)
}