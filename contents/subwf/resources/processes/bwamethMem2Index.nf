process BwamethMem2Index {
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
          path("${D}"), path("${D}.amb"), path("${D}.ann"), path("${D}.bwt.2bit.64"),
          path("${D}.pac"), path("${D}.0123")

    shell:
    alignerIndexDir = "bwameth/index"
    S = "${genomeReference}.bwameth.c2t"
    D = "${alignerIndexDir}/${prefix}.c2t"

"""
mkdir -p ${alignerIndexDir}
bwameth.py index-mem2 ${genomeReference}
mv ${S} ${D}
mv ${S}.amb ${D}.amb
mv ${S}.ann ${D}.ann
mv ${S}.bwt.2bit.64 ${D}.bwt.2bit.64
mv ${S}.pac ${D}.pac
mv ${S}.0123 ${D}.0123
"""

    stub:
    alignerIndexDir = "bwameth/index"
    D = "${alignerIndexDir}/${prefix}.c2t"
"""
mkdir -p ${alignerIndexDir}
touch ${D} ${D}.amb ${D}.ann ${D}.bwt.2bit.64 ${D}.pac ${D}.0123
"""
}