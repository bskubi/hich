def buildCmd(id, pairs, singleCell, maxMismatch, method, pairtoolsDedupParams, cpus) {
    if (!pairtoolsDedupParams) pairtoolsDedupParams = []
    if (singleCell) pairtoolsDedupParams += ["--extra-col-pair cellID cellID --backend cython --chunksize 1000000"]
    if (maxMismatch != null) pairtoolsDedupParams += ["--max-mismatch ${maxMismatch}"]
    if (method != null) pairtoolsDedupParams += ["--method ${method}"]
    output = "${id}.dedup.pairs.gz"

    cmd = "pairtools dedup --output '${output}' ${pairtoolsDedupParams.join(' ')}  --nproc-in ${cpus} --nproc-out ${cpus} '${pairs}'"
    logMap = [task: "DEDUP_PAIRS", output: output, input: [id: id, pairs: pairs, singleCell: singleCell, maxMismatch: maxMismatch, method: method, pairtoolsDedupParams: pairtoolsDedupParams]]

    return [cmd, logMap, output]
}