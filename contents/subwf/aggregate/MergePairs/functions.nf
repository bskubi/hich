def buildCmd(id, toMerge, cpus) {
    def merged = "${id}.merged.pairs.gz"
    toMerge = toMerge instanceof List ? toMerge : [toMerge]
    def mergeList = toMerge.collect({"'${it}'"})
    def logMap = [task: "MERGE_PAIRS", input: [id: id, toMerge: toMerge], output: [merged: merged]]
    def cmd = "pairtools merge --output '${merged}' --nproc-in ${cpus} --nproc-out ${cpus} ${mergeList.join(' ')}"
    def stubCmd = "touch '${merged}'"
    return [merged, mergeList, logMap, cmd, stubCmd]
}