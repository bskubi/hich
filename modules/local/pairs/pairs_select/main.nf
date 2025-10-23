process PAIRS_SELECT {
    container 'bskubi/hich_tools:latest'
    
    input:
    tuple val(id), path(pairs), val(config)

    output:
    tuple val(id), path(pairs), emit: pairs

    stub:
    pairs = "${id}.pairs_select.gz"
    """
    touch '${pairs}'
    """
}