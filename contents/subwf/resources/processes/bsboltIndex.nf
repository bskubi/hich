process BSBoltIndex {
    container params.general.hichContainer
    publishDir params.general.publish.bsboltIndex ? params.general.publish.bsboltIndex : "results",
               saveAs: {params.general.publish.bsboltIndex ? it : null},
               mode: params.general.publish.mode
    label 'whenLocal_allConsuming'
    label 'index'



    input:
    tuple val(genomeReferenceString), path(genomeReference), val(alignerIndexPrefix)

    output:
    tuple val(genomeReferenceString),
          path(alignerIndexDir),
          val(alignerIndexPrefix),
          path("${alignerIndexDir}/${prefix}.fa"),
          path("${alignerIndexDir}/${prefix}.fa.ann"),
          path("${alignerIndexDir}/${prefix}.fa.amb"),
          path("${alignerIndexDir}/${prefix}.fa.opac"),
          path("${alignerIndexDir}/${prefix}.fa.pac"),
          path("${alignerIndexDir}/${prefix}.fa.bwt"),
          path("${alignerIndexDir}/${prefix}.fa.sa")

    shell:
    alignerIndexDir = "bsbolt/index/${alignerIndexPrefix}"
    prefix = "BSB_ref"
    "python -m bsbolt Index -IA -G '${genomeReference}' -DB '${alignerIndexDir}'"

    stub:
    alignerIndexDir = "bsbolt/index/${alignerIndexPrefix}"
    prefix = "BSB_ref"
    "mkdir -p '${alignerIndexDir}' && cd '${alignerIndexDir}' && touch '${prefix}.fa' '${prefix}.fa.amb' '${prefix}.fa.ann' '${prefix}.fa.opac' '${prefix}.fa.pac' '${prefix}.fa.bwt' '${prefix}.fa.sa'"
}