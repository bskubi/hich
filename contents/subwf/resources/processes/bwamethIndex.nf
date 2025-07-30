process BwamethIndex {
    publishDir params.general.publish.bwamethIndex ? params.general.publish.bwamethIndex : "results",
               saveAs: {params.general.publish.bwamethIndex ? it : null},
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
          path("${D}"), path("${D}.amb"), path("${D}.ann"), path("${D}.bwt"),
          path("${D}.pac"), path("${D}.sa")

    shell:
    // TODO: Add withLog/stubLog pattern
    alignerIndexDir = "bwameth/index"
    S = "${genomeReference}.bwameth.c2t"
    D = "${alignerIndexDir}/${prefix}.c2t"

"""
mkdir -p ${alignerIndexDir}
bwameth.py index ${genomeReference}
mv ${S} ${D}
mv ${S}.amb ${D}.amb
mv ${S}.ann ${D}.ann
mv ${S}.bwt ${D}.bwt
mv ${S}.pac ${D}.pac
mv ${S}.sa ${D}.sa
"""

    stub:
    alignerIndexDir = "bwameth/index"
    D = "${alignerIndexDir}/${prefix}.c2t"
"""
mkdir -p ${alignerIndexDir}
touch ${D} ${D}.amb ${D}.ann ${D}.bwt ${D}.pac ${D}.sa
"""
}