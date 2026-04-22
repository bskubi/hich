process EXTRACT_3D_PAIRS {
    publishDir "results/3d_pairs/", mode: "link"
    tag "$ID"

    input:
    tuple val(ID), path(sam), path(config), path(chromsizes), val(addColumns)

    output:
    tuple val(ID), path(observations), path(statistics)

    script:
    observations = "obs/${ID}"
    statistics = "stats/${ID}"
    addColumns = addColumns ? "--add-columns ${addColumns}" : ""
    """
    hich extract 3d-pairs \
        --stats-dir '${statistics}' '${config}' '${sam}' '${observations}' \
        "pairtools parse2 -c '${chromsizes}' --nproc-in ${task.cpus} --nproc-out ${task.cpus} \
        ${addColumns} --flip --drop-seq --drop-sam"
    """

    stub:
    observations = "obs/${ID}"
    statistics = "stats/${ID}"
    """
    mkdir -p ${observations}
    mkdir -p ${statistics}
    """
}