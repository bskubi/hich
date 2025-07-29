process BwaMemIndex {
    publishDir params.general.publish.bwaIndex ? params.general.publish.bwaIndex : "results",
               saveAs: {params.general.publish.bwaIndex ? it : null},
               mode: params.general.publish.mode
    label 'whenLocal_allConsuming'
    label 'index'
    memory {15.GB + 10.GB * (task.attempt - 1)}
    conda "$projectDir/env/dev_env.yml"
    container params.general.alignmentContainer


    input:
    tuple val(genomeReferenceString), path(genomeReference), val(prefix)

    output:
    tuple val(genomeReferenceString), path(alignerIndexDir), val(prefix), 
          path("${alignerIndexDir}/${prefix}.ann"), path("${alignerIndexDir}/${prefix}.amb"),
          path("${alignerIndexDir}/${prefix}.pac"), path("${alignerIndexDir}/${prefix}.bwt"), 
          path("${alignerIndexDir}/${prefix}.sa")

    shell:
    alignerIndexDir = "bwa/index"
    "bwa index -p '${prefix}' ${genomeReference} && mkdir -p '${alignerIndexDir}' && mv '${prefix}.amb' '${prefix}.ann' '${prefix}.pac' '${prefix}.bwt' '${prefix}.sa' '${alignerIndexDir}'"

    stub:
    alignerIndexDir = "bwa/index"
    "mkdir -p '${alignerIndexDir}' && cd '${alignerIndexDir}' && touch '${prefix}.amb' '${prefix}.ann' '${prefix}.pac' '${prefix}.bwt' '${prefix}.sa'"
}