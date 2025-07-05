include {emptyOnLastStep; skip} from '../util/cli.nf'
include {keyUpdate} from '../util/keyUpdate.nf'
include {isExistingFile} from '../util/files.nf'
include {withLog; stubLog} from '../util/logs.nf'

process BwaMem2Index {
    publishDir params.general.publish.bwaMem2Index ? params.general.publish.bwaMem2Index : "results",
               saveAs: {params.general.publish.bwaMem2Index ? it : null},
               mode: params.general.publish.mode
    label 'whenLocal_allConsuming'
    label 'index'

    input:
    tuple val(genomeReferenceString), path(genomeReference), val(prefix)

    output:
    tuple val(genomeReferenceString), path("bwa-mem2/index"), val(prefix), path("bwa-mem2/index/${prefix}.0123"), path("bwa-mem2/index/${prefix}.amb"),
          path("bwa-mem2/index/${prefix}.ann"), path("bwa-mem2/index/${prefix}.bwt.2bit.64"),
          path("bwa-mem2/index/${prefix}.pac")

    shell:
    alignerIndexDir = "bwa-mem2/index"
    "bwa-mem2 index -p '${prefix}' '${genomeReference}' && mkdir -p '${alignerIndexDir}' && mv '${prefix}.0123' '${prefix}.amb' '${prefix}.ann' '${prefix}.bwt.2bit.64' '${prefix}.pac' '${alignerIndexDir}'"

    stub:
    alignerIndexDir = "bwa-mem2/index"
    "mkdir -p '${alignerIndexDir}' && cd '${alignerIndexDir}' && touch '${prefix}.0123' '${prefix}.amb' '${prefix}.ann' '${prefix}.bwt.2bit.64' '${prefix}.pac'"
}

process BwaMemIndex {
    publishDir params.general.publish.bwaIndex ? params.general.publish.bwaIndex : "results",
               saveAs: {params.general.publish.bwaIndex ? it : null},
               mode: params.general.publish.mode
    label 'whenLocal_allConsuming'
    label 'index'
    memory {20.GB + 10.GB * (task.attempt - 1)}
    


    input:
    tuple val(genomeReferenceString), path(genomeReference), val(prefix)

    output:
    tuple val(genomeReferenceString), path(alignerIndexDir), val(prefix), path("${alignerIndexDir}/${prefix}.ann"), path("${alignerIndexDir}/${prefix}.amb"),
          path("${alignerIndexDir}/${prefix}.pac"), path("${alignerIndexDir}/${prefix}.bwt"), path("${alignerIndexDir}/${prefix}.sa")

    shell:
    alignerIndexDir = "bwa/index"
    "bwa index -p '${prefix}' ${genomeReference} && mkdir -p '${alignerIndexDir}' && mv '${prefix}.amb' '${prefix}.ann' '${prefix}.pac' '${prefix}.bwt' '${prefix}.sa' '${alignerIndexDir}'"

    stub:
    alignerIndexDir = "bwa/index"
    "mkdir -p '${alignerIndexDir}' && cd '${alignerIndexDir}' && touch '${prefix}.amb' '${prefix}.ann' '${prefix}.pac' '${prefix}.bwt' '${prefix}.sa'"
}

process BSBoltIndex {
    container params.general.hichContainer
    publishDir params.general.publish.bsboltIndex ? params.general.publish.bsboltIndex : "results",
               saveAs: {params.general.publish.bsboltIndex ? it : null},
               mode: params.general.publish.mode
    label 'whenLocal_allConsuming'
    label 'index'
    memory {20.GB + 10.GB * (task.attempt - 1)}
    


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


workflow AlignerIndex {

    take:
        samples
    
    main:
    samples
        | filter {!skip("alignerIndex") && it.datatype == "fastq" && it.aligner == "bwa-mem2"}
        | filter {!isExistingFile(it.alignerIndexDir)}
        | map{it.alignerIndexPrefix = it.alignerIndexPrefix ?: it.assembly; it}
        | map{tuple(it.genomeReference, it.genomeReference, it.alignerIndexPrefix)}
        | unique
        | BwaMem2Index
        | map{genomeReference, alignerIndexDir, alignerIndexPrefix, prefix_0123, prefix_amb, prefix_ann, prefix_bwt_2bit_64, prefix_pac ->
                [genomeReference: file(genomeReference), alignerIndexDir: alignerIndexDir, alignerIndexPrefix: alignerIndexPrefix]}
        | set{resultBwaMem2Index}

    samples
        | filter {!skip("alignerIndex") && it.datatype == "fastq" && it.aligner == "bwa"}
        | filter {!isExistingFile(it.alignerIndexDir)}
        | map{it.alignerIndexPrefix = it.alignerIndexPrefix ?: it.assembly; it}
        | map{tuple(it.genomeReference, it.genomeReference, it.alignerIndexPrefix)}
        | unique
        | BwaMemIndex
        | map{genomeReference, alignerIndexDir, alignerIndexPrefix, prefix_ann, prefix_amb, prefix_pac, prefix_bwt, prefix_sa ->
                [genomeReference: file(genomeReference), alignerIndexDir: alignerIndexDir, alignerIndexPrefix: alignerIndexPrefix]}
        | set{resultBwaMemIndex}

    samples
        | filter {!skip("alignerIndex") && it.datatype == "fastq" && it.aligner == "bsbolt"}
        | filter {!isExistingFile(it.alignerIndexDir)}
        | map{it.alignerIndexPrefix = it.alignerIndexPrefix ?: it.assembly; it}
        | map{tuple(it.genomeReference, it.genomeReference, it.alignerIndexPrefix)}
        | unique
        | BSBoltIndex
        | map{genomeReference, alignerIndexDir, alignerIndexPrefix,
              prefix_fa, prefix_ann, prefix_amb, prefix_opac, prefix_pac, prefix_bwt, prefix_sa ->
                [genomeReference: file(genomeReference), alignerIndexDir: alignerIndexDir, alignerIndexPrefix: alignerIndexPrefix]}
        | set{resultBSBoltIndex}

    keyUpdate(samples, resultBwaMem2Index, "genomeReference") | set{samples}
    keyUpdate(samples, resultBwaMemIndex, "genomeReference") | set{samples}
    keyUpdate(samples, resultBSBoltIndex, "genomeReference") | set{samples}

    samples = emptyOnLastStep("alignerIndex", samples)

    emit:
        samples

}

