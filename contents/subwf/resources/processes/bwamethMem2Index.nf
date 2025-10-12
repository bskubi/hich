process BWAMETH_MEM2_INDEX {
    publishDir params.general.publish.alignerIndex ?: "results",
               saveAs: {params.general.publish.alignerIndex ? it : null},
               mode: params.general.publish.mode
    label 'fullNodeLongTime'
    memory {15.GB + 10.GB * (task.attempt - 1)}
    conda "$projectDir/env/dev_env.yml"
    container params.general.alignmentContainer

    input:
    tuple val(genomeReferenceString), path(genomeReference), val(prefix)

    output:
    tuple val(genomeReferenceString), path(alignerIndexDir), val(prefix), 
          path(iC2t), path(iAmb), path(iAnn), path(iBwt2Bit64),
          path(iPac), path(i0123)

    shell:
    alignerIndexDir = "bwameth-mem2"
    S = "${genomeReference}.bwameth.c2t"
    D = "${alignerIndexDir}/${prefix}.c2t"
    iC2t = D
    iAmb = "${D}.amb"
    iAnn = "${D}.ann"
    iBwt2Bit64 = "${D}.bwt.2bit.64"
    iPac = "${D}.pac"
    i0123 = "${D}.0123"
    """
    mkdir -p '${alignerIndexDir}'
    bwameth.py index-mem2 '${genomeReference}'
    mv '${S}' '${iC2t}'
    mv '${S}.amb' '${iAmb}'
    mv '${S}.ann' '${iAnn}'
    mv '${S}.bwt.2bit.64' '${iBwt2Bit64}'
    mv '${S}.pac' '${iPac}'
    mv '${S}.0123' '${i0123}'
    """

    stub:
    alignerIndexDir = "bwameth-mem2"
    S = "${genomeReference}.bwameth.c2t"
    D = "${alignerIndexDir}/${prefix}.c2t"
    iC2t = D
    iAmb = "${D}.amb"
    iAnn = "${D}.ann"
    iBwt2Bit64 = "${D}.bwt.2bit.64"
    iPac = "${D}.pac"
    i0123 = "${D}.0123"
    """
    mkdir -p '${alignerIndexDir}'
    touch '${iC2t}' '${iAmb}' '${iAnn}' '${iBwt2Bit64}' '${iPac}' '${i0123}'
    """
}