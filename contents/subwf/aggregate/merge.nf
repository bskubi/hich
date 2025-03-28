include {withLog; stubLog} from '../util/logs.nf'

process Merge {
    label 'pairs'
    
    input:
    tuple val(id), path(toMerge)

    output:
    tuple val(id), path(merged)

    shell:
    merged = "${id}.merged.pairs.gz"
    toMerge = toMerge.collect { "'${it}'" }
    cmd = "pairtools merge --output '${merged}'  --nproc-in ${task.cpus} --nproc-out ${task.cpus} ${toMerge.join(' ')}"
    logMap = [task: "PairtoolsMerge", input: [id: id, toMerge: toMerge], output: [merged: merged]]
    withLog(cmd, logMap)

    stub:
    merged = "${id}.merged.pairs.gz"
    stub = "touch '${merged}'"
    toMerge = toMerge.collect { "'${it}'" }
    cmd = "pairtools merge --output '${merged}'  --nproc-in ${task.cpus} --nproc-out ${task.cpus} ${toMerge.join(' ')}"
    logMap = [task: "PairtoolsMerge", input: [id: id, toMerge: toMerge], output: [merged: merged]]
    stubLog(stub, cmd, logMap)  
}