process PAIRS_DEDUP {
    container 'bskubi/hich_tools:latest'
    
    input:
    tuple val(id), path(pairs), val(config)

    output:
    tuple val(id), path(pairs), emit: pairs

    stub:
    pairs = "${id}.pairs_dedup.gz"
    """
    touch '${pairs}'
    """
}