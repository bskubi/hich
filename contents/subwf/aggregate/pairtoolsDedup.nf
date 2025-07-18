include {withLog; stubLog} from '../util/logs.nf'

process PairtoolsDedup {
    publishDir params.general.publish.parse ? params.general.publish.dedup : "results",
               saveAs: {params.general.publish.dedup ? it : null},
               mode: params.general.publish.mode

    label 'pairs'
    conda "$projectDir/env/dev_env.yml"
    container params.general.hichContainer

    input:
    tuple val(id), path(pairs), val(singleCell), val(maxMismatch), val(method), val(pairtoolsDedupParams)

    output:
    tuple val(id), path(deduplicated)

    shell:
    if (!pairtoolsDedupParams) pairtoolsDedupParams = []
    if (singleCell) pairtoolsDedupParams += ["--extra-col-pair cellID cellID --backend cython --chunksize 1000000"]
    if (maxMismatch != null) pairtoolsDedupParams += ["--max-mismatch ${maxMismatch}"]
    if (method != null) pairtoolsDedupParams += ["--method ${method}"]
    deduplicated = "${id}.dedup.pairs.gz"

    cmd = "pairtools dedup --output '${deduplicated}' ${pairtoolsDedupParams.join(' ')}  --nproc-in ${task.cpus} --nproc-out ${task.cpus} '${pairs}'"
    logMap = [task: "PairtoolsDedup", output: deduplicated, input: [id: id, pairs: pairs, singleCell: singleCell, maxMismatch: maxMismatch, method: method, pairtoolsDedupParams: pairtoolsDedupParams]]
    withLog(cmd, logMap)

    stub:
    deduplicated = "${id}.dedup.pairs.gz"
    stub = "touch '${deduplicated}'"
    if (!pairtoolsDedupParams) pairtoolsDedupParams = []
    if (singleCell) pairtoolsDedupParams += ["--extra-col-pair cellID cellID --backend cython --chunksize 1000000"]
    if (maxMismatch != null) pairtoolsDedupParams += ["--max-mismatch ${maxMismatch}"]
    if (method != null) pairtoolsDedupParams += ["--method ${method}"]
    deduplicated = "${id}.dedup.pairs.gz"

    cmd = "pairtools dedup --output '${deduplicated}' ${pairtoolsDedupParams.join(' ')}  --nproc-in ${task.cpus} --nproc-out ${task.cpus} '${pairs}'"
    logMap = [task: "PairtoolsDedup", output: deduplicated, input: [id: id, pairs: pairs, singleCell: singleCell, maxMismatch: maxMismatch, method: method, pairtoolsDedupParams: pairtoolsDedupParams]]
    stubLog(stub, cmd, logMap)
}
