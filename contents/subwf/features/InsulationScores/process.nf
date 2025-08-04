include {withLog; stubLog} from '../../util/logs.nf'
include {buildCmd} from './functions.nf'

process INSULATION_SCORES {
    publishDir "results/insulation",
               mode: params.general.publish.mode

    conda "$projectDir/env/dev_env.yml"
    container params.general.hichContainer
    label 'features'
    tag "$id"

    input:
    tuple val(id), path(mcool), val(resolution), val(cooltoolsInsulationParams), val(window)

    output:
    tuple val(id), path(tsv), path(bw)

    shell:
    (cmd, logMap, tsv, bw) = buildCmd(id, mcool, resolution, cooltoolsInsulationParams, window)
    withLog(cmd, logMap)

    stub:
    (cmd, logMap, tsv, bw) = buildCmd(id, mcool, resolution, cooltoolsInsulationParams, window)
    stub = "touch '${tsv}' '${bw}'"
    stubLog(stub, cmd, logMap)
}