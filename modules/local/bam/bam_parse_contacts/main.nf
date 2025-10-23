process BAM_PARSE_CONTACTS {
    container 'bskubi/hich_tools:latest'
    
    input:
    tuple val(id), path(bam), val(config_bam_parse_contacts)

    output:
    tuple val(id), path(pairs), emit: pairs

    shell:
    pairs = "${id}.bam_parse_contacts.pairs.gz"
    """
    touch '${pairs}'
    """
}