process MergeTechrepsToBioreps {
    container params.general.hichContainer
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


workflow Merge {
    take:
    samples

    main:

    samples
    | filter{it.aggregateLevel == "techrep"}
    
    // Merge the pairs files.
    groupHashMap(merge.YES, levelParams.mergeGroupIdentifiers)
        | map{coalesce(columns(it, ["dropNull":true]))}
        | map{tuple(constructIdentifier(coalesce(it, "_drop")), it.latestPairs)}
        | MergeTechrepsToBioreps
        | map{[id:it[0], pairs:it[1], latest:it[1], latestPairs:it[1]]}
        | set{fromMerge}

    samples = emptyOnLastStep("merge", samples)
    emit:
    samples
}