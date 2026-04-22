process MERGE_DEDUP {
    publishDir = "results/merge_dedup/"
    tag "$ID"

    input:
    tuple val(ID), path(sams), val(dedupTag), path(configEditSAM)

    output:
    tuple val(ID), path(samDedup)

    script:
    statsOutput = "${ID}.stats.json"
    samDedup = "${ID}.dedup.bam"
    sams = [sams].flatten().collect{it.toString()}.join(" ")

    nProcsDedup = "--n-procs ${task.cpus}"
    dedupTag = dedupTag ? "--barcode-tag ${dedupTag}" : ""
    configEditSAM = configEditSAM ? "--config-edit-sam '${configEditSAM}'" : ""
    samOutput = "--sam-output ${samDedup}"

    """
    samtools cat ${sams} | dedup.py ${nProcsDedup} ${dedupTag} ${configEditSAM} ${samOutput}
    """

    stub:
    samDedup = "${ID}.dedup.bam"
    sams = [sams].flatten().collect{it.toString()}.join(" ")
    """
    touch ${samDedup}
    """
}