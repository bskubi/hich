include {withLog; stubLog} from '../../util/logs.nf'
include {buildCmd} from './functions.nf'


process HICREP_COMBINATIONS {
    publishDir "results/hicrep", mode: params.general.publish.mode
    tag "$planName"

    conda "$projectDir/env/dev_env.yml"
    container params.general.hichContainer
    
    input:
    tuple val(planName), path(mcools), val(hicrep_combinations_opts)

    output:
    path(output)

    shell:
    (cmd, logMap, output) = buildCmd(planName, mcools, hicrep_combinations_opts, task.cpus)
    withLog(cmd, logMap)

    stub:
    (cmd, logMap, output) = buildCmd(planName, mcools, hicrep_combinations_opts, task.cpus)
    stub = "touch '${output}'"
    stubLog(stub, cmd, logMap)
}