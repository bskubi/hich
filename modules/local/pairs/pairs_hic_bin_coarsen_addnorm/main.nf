process PAIRS_HIC_BIN_COARSEN_ADDNORM {
    container 'bskubi/hich_juicer:latest'

    input:
    tuple val(id), path(pairs), val(config)

    output:
    tuple val(id), path(hic), emit: hic

    stub:
    hic = "${id}.hic"
    """
    touch '${hic}'
    """
}