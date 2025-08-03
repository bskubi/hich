include {withLog; stubLog} from '../../util/logs.nf'
include {buildCmd} from './functions.nf'

process MERGE_PAIRS {
    label 'pairs'
    conda "$projectDir/env/dev_env.yml"
    container params.general.hichContainer
    
    input:
    tuple val(id), path(toMerge)

    output:
    tuple val(id), path(merged)

    shell:
    (merged, mergeList, logMap, cmd, stubCmd) = buildCmd(id, toMerge, task.cpus)
    withLog(cmd, logMap)

    stub:
    (merged, mergeList, logMap, cmd, stubCmd) = buildCmd(id, toMerge, task.cpus)
    withLog(cmd, logMap, stubCmd)
}
