process PAIRS_LABEL {
    container 'bskubi/hich_tools:latest'
    
    input:
    tuple val(id), path(pairs), val(config)

    output:
    tuple val(id), path(pairs), emit: pairs

    stub:
    pairs = "${id}.pairs_label.gz"
    """
    touch '${pairs}'
    """
}