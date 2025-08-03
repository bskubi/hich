include {withLog; stubLog} from '../../util/logs.nf'
include {buildCmdMcoolToHic; buildCmdHicToMcool} from './functions.nf'
process MCOOL_TO_HIC {
    publishDir params.general.publish.hic ? params.general.publish.hic : "results",
               saveAs: {params.general.publish.hic ? it : null},
               mode: params.general.publish.mode
    
    label 'convertMcoolToHic'
    conda "$projectDir/env/dev_env.yml"
    container params.general.hictkContainer

    input:
    tuple val(id), path(mcoolFile)

    output:
    tuple val(id), path(output)

    shell:
    (cmd, logMap, output) = buildCmdMcoolToHic(id, mcoolFile, task.cpus)
    withLog(cmd, logMap)

    stub:
    (cmd, logMap, output) = buildCmdMcoolToHic(id, mcoolFile, task.cpus) 
    stub = "touch '${output}'"
    stubLog(stub, cmd, logMap)

}

process HIC_TO_MCOOL {
    publishDir params.general.publish.mcool ? params.general.publish.mcool : "results",
               saveAs: {params.general.publish.mcool ? it : null},
               mode: params.general.publish.mode
    
    label 'convertHicToMcool'
    cpus 2
    maxForks 1  // Necessary due to unexplained bug in hictk
    conda "$projectDir/env/dev_env.yml"
    container params.general.hictkContainer

    input:
    tuple val(id), path(hicFile)

    output:
    tuple val(id), path(output)

    shell:
    (cmd, logMap, output) = buildCmdHicToMcool(id, hicFile)
    withLog(cmd, logMap)

    stub:
    (cmd, logMap, output) = buildCmdHicToMcool(id, hicFile)
    stub = "touch '${output}'"
    stubLog(stub, cmd, logMap)

}