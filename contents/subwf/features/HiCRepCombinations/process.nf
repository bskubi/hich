include {withLog; stubLog} from '../../util/logs.nf'
include {buildCmd} from './functions.nf'


process HICREP_COMBINATIONS {
    publishDir "results/hicrep", mode: params.general.publish.mode
    tag "$planName"

    conda "$projectDir/env/dev_env.yml"
    container params.general.hichContainer
    
    input:
    tuple val(planName), path(mcools), val(resolutions), val(chroms), val(exclude), val(chromFilter), val(h), val(dBPMax), val(bDownSample)

    output:
    path(output)

    shell:
    (cmd, logMap, output) = buildCmd(planName, mcools, resolutions, chroms, exclude, chromFilter, h, dBPMax, bDownSample)
    withLog(cmd, logMap)

    stub:
    (cmd, logMap, output) = buildCmd(planName, mcools, resolutions, chroms, exclude, chromFilter, h, dBPMax, bDownSample)
    stub = "touch '${output}'"
    stubLog(stub, cmd, logMap)
}