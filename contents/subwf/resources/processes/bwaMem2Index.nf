process BwaMem2Index {
    publishDir params.general.publish.bwaMem2Index ? params.general.publish.bwaMem2Index : "results",
               saveAs: {params.general.publish.bwaMem2Index ? it : null},
               mode: params.general.publish.mode
    label 'whenLocal_allConsuming'
    label 'index'
    conda "$projectDir/env/dev_env.yml"
    container params.general.alignmentContainer

    input:
    tuple val(genomeReferenceString), path(genomeReference), val(prefix)

    output:
    tuple val(genomeReferenceString), path("bwa-mem2/index"), val(prefix), 
          path("bwa-mem2/index/${prefix}.0123"), path("bwa-mem2/index/${prefix}.amb"),
          path("bwa-mem2/index/${prefix}.ann"), path("bwa-mem2/index/${prefix}.bwt.2bit.64"),
          path("bwa-mem2/index/${prefix}.pac")

    shell:
    alignerIndexDir = "bwa-mem2/index"
    "bwa-mem2 index -p '${prefix}' '${genomeReference}' && mkdir -p '${alignerIndexDir}' && mv '${prefix}.0123' '${prefix}.amb' '${prefix}.ann' '${prefix}.bwt.2bit.64' '${prefix}.pac' '${alignerIndexDir}'"

    stub:
    alignerIndexDir = "bwa-mem2/index"
    "mkdir -p '${alignerIndexDir}' && cd '${alignerIndexDir}' && touch '${prefix}.0123' '${prefix}.amb' '${prefix}.ann' '${prefix}.bwt.2bit.64' '${prefix}.pac'"
}