process COOL_COARSEN_ADDNORM {
    container 'bskubi/hich_cooler:latest'

    input:
    tuple val(id), path(cool), val(config)

    output:
    tuple val(id), path(mcool), emit: mcool

    stub:
    mcool = "${id}.mcool"
    """
    touch '${mcool}'
    """
}