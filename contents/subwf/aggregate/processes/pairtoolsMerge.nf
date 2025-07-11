include {withLog; stubLog} from '../../util/logs.nf'

def buildCmd(id, toMerge, cpus) {
    def merged = "${id}.merged.pairs.gz"
    toMerge = toMerge instanceof List ? toMerge : [toMerge]
    def mergeList = toMerge.collect({"'${it}'"})
    def logMap = [task: "Merge", input: [id: id, toMerge: toMerge], output: [merged: merged]]
    def cmd = "pairtools merge --output '${merged}' --nproc-in ${cpus} --nproc-out ${cpus} ${mergeList.join(' ')}"
    def stubCmd = "touch '${merged}'"
    return [merged, mergeList, logMap, cmd, stubCmd]
}

process PairtoolsMerge {
    label 'pairs'
    conda "$projectDir/env/dev_env.yml"
    
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
