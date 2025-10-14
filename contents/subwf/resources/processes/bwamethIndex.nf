process BWAMETH_INDEX {
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
          path(iC2t), path(iAmb), path(iAnn), path(iBwt),
          path(iPac), path(iSa)

    shell:
    alignerIndexDir = "bwameth"
    S = "${genomeReference}.bwameth.c2t"
    D = "${alignerIndexDir}/${prefix}.c2t"
    iC2t = D
    iAmb = "${D}.amb"
    iAnn = "${D}.ann"
    iBwt = "${D}.bwt"
    iPac = "${D}.pac"
    iSa = "${D}.sa"

    """
    mkdir -p '${alignerIndexDir}'
    bwameth.py index '${genomeReference}'
    mv ${S} ${iC2t}
    mv ${S}.amb ${iAmb}
    mv ${S}.ann ${iAnn}
    mv ${S}.bwt ${iBwt}
    mv ${S}.pac ${iPac}
    mv ${S}.sa ${iSa}
    """

    stub:
    alignerIndexDir = "bwameth"
    S = "${genomeReference}.bwameth.c2t"
    D = "${alignerIndexDir}/${prefix}.c2t"
    iC2t = D
    iAmb = "${D}.amb"
    iAnn = "${D}.ann"
    iBwt = "${D}.bwt"
    iPac = "${D}.pac"
    iSa = "${D}.sa"
    """
    mkdir -p '${alignerIndexDir}'
    touch '${iC2t}' '${iAmb}' '${iAnn}' '${iBwt}' '${iPac}' '${iSa}'
    """
}