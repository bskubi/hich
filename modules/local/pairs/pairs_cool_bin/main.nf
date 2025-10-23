process PAIRS_COOL_BIN {
    container 'bskubi/hich_cooler:latest'

    input:
    tuple val(id), path(pairs), val(config)

    output:
    tuple val(id), path(cool), emit: cool

    stub:
    cool = "${id}.cool"
    """
    touch '${cool}'
    """
}