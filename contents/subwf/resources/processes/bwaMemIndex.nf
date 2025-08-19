process BwaMemIndex {
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
          path("${alignerIndexDir}/${iAnn}"), path("${alignerIndexDir}/${iAmb}"),
          path("${alignerIndexDir}/${iPac}"), path("${alignerIndexDir}/${iBwt}"), 
          path("${alignerIndexDir}/${iSa}")

    shell:
    alignerIndexDir = "bwa"
    iAnn = "${prefix}.ann"
    iAmb = "${prefix}.amb"
    iPac = "${prefix}.pac"
    iBwt = "${prefix}.bwt"
    iSa = "${prefix}.sa"
    """
    mkdir -p '${alignerIndexDir}'
    bwa index -p '${prefix}' '${genomeReference}'
    mv '${iAnn}' '${iAmb}' '${iPac}' '${iBwt}' '${iSa}' '${alignerIndexDir}'
    """

    stub:
    alignerIndexDir = "bwa"
    iAnn = "${prefix}.ann"
    iAmb = "${prefix}.amb"
    iPac = "${prefix}.pac"
    iBwt = "${prefix}.bwt"
    iSa = "${prefix}.sa"
    """
    mkdir -p '${alignerIndexDir}'
    cd '${alignerIndexDir}'
    touch '${iAnn}' '${iAmb}' '${iPac}' '${iBwt}' '${iSa}'
    """
}