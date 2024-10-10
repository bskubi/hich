include {emptyOnLastStep; isExistingFile; pack; skip} from './extraops.nf'

process BwaMem2Index {
    conda "bioconda::bwa-mem2"
    container "bskubi/hich:latest"
    publishDir params.general.publish.bwa_mem2Index ? params.general.publish.bwa_mem2Index : "results",
               saveAs: {params.general.publish.bwa_mem2Index ? it : null},
               mode: params.general.publish.mode
    label 'whenLocal_allConsuming'
    label 'index'

    input:
    tuple path(genomeReference), val(prefix)

    output:
    tuple val(genomeReference), path("bwa-mem2/index"), val(prefix), path("bwa-mem2/index/${prefix}.0123"), path("bwa-mem2/index/${prefix}.amb"),
          path("bwa-mem2/index/${prefix}.ann"), path("bwa-mem2/index/${prefix}.bwt.2bit.64"),
          path("bwa-mem2/index/${prefix}.pac")

    shell:
    alignerIndexDir = "bwa-mem2/index"
    "bwa-mem2 index -p ${prefix} ${genomeReference} && mkdir -p ${alignerIndexDir} && mv ${prefix}.0123 ${prefix}.amb ${prefix}.ann ${prefix}.bwt.2bit.64 ${prefix}.pac ${alignerIndexDir}"

    stub:
    alignerIndexDir = "bwa-mem2/index"
    "mkdir -p ${alignerIndexDir} && cd ${alignerIndexDir} && touch ${prefix}.0123 ${prefix}.amb ${prefix}.ann ${prefix}.bwt.2bit.64 ${prefix}.pac"
}

process BwaMemIndex {
    conda "bioconda::bwa"
    container "bskubi/hich:latest"
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
    "bwa index -p ${prefix} ${genomeReference} && mkdir -p ${alignerIndexDir} && mv ${prefix}.amb ${prefix}.ann ${prefix}.pac ${prefix}.bwt ${prefix}.sa ${alignerIndexDir}"

    stub:
    alignerIndexDir = "bwa/index"
    "mkdir -p ${alignerIndexDir} && cd ${alignerIndexDir} && touch ${prefix}.amb ${prefix}.ann ${prefix}.pac ${prefix}.bwt ${prefix}.sa"
}

process BSBoltIndex {
    container "bskubi/hich:latest"
    publishDir params.general.publish.indexBSBolt ? params.general.publish.indexBSBolt : "results",
               saveAs: {params.general.publish.indexBSBolt ? it : null},
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
    alignerIndexDir = "bsbolt/index"
    "bsbolt Index -G ${genomeReference} && mkdir -p ${alignerIndexDir} && mv ${prefix}.amb ${prefix}.ann ${prefix}.pac ${prefix}.bwt ${prefix}.sa ${alignerIndexDir}"

    stub:
    alignerIndexDir = "bsbolt/index"
    "mkdir -p ${alignerIndexDir} && cd ${alignerIndexDir} && touch ${prefix}.amb ${prefix}.ann ${prefix}.pac ${prefix}.bwt ${prefix}.sa"
}


workflow AlignerIndex {

    take:
        samples
    
    main:
    samples
        | filter {!skip("alignerIndex") && it.datatype == "fastq" && it.aligner == "bwa-mem2"}
        | filter {!isExistingFile(it.alignerIndexDir)}
        | map{it.alignerIndexPrefix = it.alignerIndexPrefix ?: it.assembly; it}
        | map{tuple(it.genomeReference, it.alignerIndexPrefix)}
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
    
    pack(samples, resultBwaMem2Index, "genomeReference") | set{samples}
    pack(samples, resultBwaMemIndex, "genomeReference") | set{samples}

    samples = emptyOnLastStep("alignerIndex", samples)

    emit:
        samples

}