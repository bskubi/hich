process PAIRS_MERGE {
    container 'bskubi/hich_tools:latest'
    
    input:
    tuple val(id), path(pairs), val(config), val(label_suffix)

    output:
    tuple val(id), path(pairs), emit: pairs

    stub:
    pairs = "${id}.pairs_merge${label_suffix}.gz"
    """
    touch '${pairs}'
    """
}