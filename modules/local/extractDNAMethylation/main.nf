process EXTRACT_DNA_METHYLATION {
    publishDir = "results/dna_methylation/"
    tag "$ID"
    
    input:
    tuple val(ID), path(sam), path(genomeReference), path(config)

    output:
    tuple val(ID), path(observations), path(statistics)

    script:
    observations = "obs/${ID}"
    statistics = "stats/${ID}"
    """
    hich extract dna-methylation \
        --limit-chroms chr1 \
        --stats-dir ${statistics} \
        --n-workers ${task.cpus} \
        --worker-n-procs 1 \
        ${config} ${sam} ${genomeReference} ${observations}
    """

    stub:
    observations = "obs/${ID}"
    statistics = "stats/${ID}"
    """
    mkdir -p ${observations}
    mkdir -p ${statistics}
    """
}