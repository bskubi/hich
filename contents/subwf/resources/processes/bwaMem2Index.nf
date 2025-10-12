process BWA_MEM2_INDEX {
    publishDir params.general.publish.alignerIndex ?: "results",
               saveAs: {params.general.publish.alignerIndex ? it : null},
               mode: params.general.publish.mode
    label 'fullNodeLongTime'
    conda "$projectDir/env/dev_env.yml"
    container params.general.alignmentContainer

    input:
    tuple val(genomeReferenceString), path(genomeReference), val(prefix)

    output:
    tuple val(genomeReferenceString), path(alignerIndexDir), val(prefix), 
          path("${alignerIndexDir}/${i0123}"), path("${alignerIndexDir}/${iAmb}"),
          path("${alignerIndexDir}/${iAnn}"), path("${alignerIndexDir}/${iBwt2Bit64}"),
          path("${alignerIndexDir}/${iPac}")

    shell:
    alignerIndexDir = "bwa-mem2"
    i0123 = "${prefix}.0123"
    iAmb = "${prefix}.amb"
    iAnn = "${prefix}.ann"
    iBwt2Bit64 = "${prefix}.bwt.2bit.64"
    iPac = "${prefix}.pac"
    """
    mkdir -p '${alignerIndexDir}'
    bwa-mem2 index -p '${prefix}' '${genomeReference}'
    mv '${i0123}' '${iAmb}' '${iAnn}' '${iBwt2Bit64}' '${iPac}' '${alignerIndexDir}'
    """

    stub:
    alignerIndexDir = "bwa-mem2"
    i0123 = "${prefix}.0123"
    iAmb = "${prefix}.amb"
    iAnn = "${prefix}.ann"
    iBwt2Bit64 = "${prefix}.bwt.2bit.64"
    iPac = "${prefix}.pac"
    """
    mkdir -p '${alignerIndexDir}'
    cd '${alignerIndexDir}'
    touch '${i0123}' '${iAmb}' '${iAnn}' '${iBwt2Bit64}' '${iPac}'
    """
}