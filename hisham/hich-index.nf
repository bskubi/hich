process index {
    tag       "${genome}"
    container "bskubi/bwa-mem2:latest"
    publishDir "resources/index/bwa-mem2", mode: 'move', overwrite: false

    input:
    tuple val(assembly_id), path(reference_fasta)

    output:
    tuple(path("${assembly_id}.0123"), 
          path("${assembly_id}.amb"),
          path("${assembly_id}.ann"),
          path("${assembly_id}.bwt.2bit.64"), 
          path("${assembly_id}.pac"))
    
    shell:
        [
            "bwa-mem2 index -p !{assembly_id} !{reference_fasta}"
        ].join(' ')

}

workflow {
    index_ch = channel.of("hg38_noalts").combine(channel.fromPath("hg38_noalts.fasta.gz"))
    index_ch | view
    index(index_ch)
}
