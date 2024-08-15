include {sourcePrefix; emptyOnLastStep} from './extraops.nf'

process BwaMem2Index {
    conda "bioconda::bwa-mem2"
    container "bskubi/hich:latest"
    publishDir params.general.publish.bwa_mem2_index ? params.general.publish.bwa_mem2_index : "results",
               saveAs: {params.general.publish.bwa_mem2_index ? it : null},
               mode: params.general.publish.mode

    input:
    tuple path(genomeReference), val(prefix)

    output:
    tuple path("bwa-mem2/index"), val(prefix), path("bwa-mem2/index/${prefix}.0123"), path("bwa-mem2/index/${prefix}.amb"),
          path("bwa-mem2/index/${prefix}.ann"), path("bwa-mem2/index/${prefix}.bwt.2bit.64"),
          path("bwa-mem2/index/${prefix}.pac")

    shell:
    "bwa-mem2 index -p ${prefix} ${genomeReference} && mkdir -p bwa-mem2/index && mv ${prefix}.0123 ${prefix}.amb ${prefix}.ann ${prefix}.bwt.2bit.64 ${prefix}.pac bwa-mem2/index"

    stub:
    "mkdir -p bwa-mem2/index && cd bwa-mem2/index && touch ${prefix}.0123 ${prefix}.amb ${prefix}.ann ${prefix}.bwt.2bit.64 ${prefix}.pac"
}

process BwaMemIndex {
    conda "bioconda::bwa"
    container "bskubi/hich:latest"
    publishDir params.general.publish.bwa_mem_index ? params.general.publish.bwa_mem_index : "results",
               saveAs: {params.general.publish.bwa_mem_index ? it : null},
               mode: params.general.publish.mode

    input:
    tuple path(genomeReference), val(prefix)

    output:
    tuple path("bwa/index"), val(prefix), path("bwa/index/${prefix}.ann"), path("bwa/index/${prefix}.amb"),
          path("bwa/index/${prefix}.pac"), path("bwa/index/${prefix}.bwt"), path("bwa/index/${prefix}.sa")

    shell:
    "bwa index -p ${prefix} ${genomeReference} && mkdir -p bwa/index && mv ${prefix}.amb ${prefix}.ann ${prefix}.pac ${prefix}.bwt ${prefix}.sa bwa/index"

    stub:
    "mkdir -p bwa/index && cd bwa/index && touch ${prefix}.amb ${prefix}.ann ${prefix}.pac ${prefix}.bwt ${prefix}.sa"
}


workflow MakeMissingIndex {

    take:
        samples
    
    main:
        sourcePrefix(
            BwaMem2Index,
            samples,
            "index_dir",
            "index_prefix",
            ["genomeReference", "index_prefix"],
            ["index_dir", "index_prefix", ".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac"],
            {["index_prefix":it.assembly]},
            "index_prefix",
            {it.datatype in ["fq", "fastq"] && it.aligner == "bwa-mem2"},
            ["keep":["index_dir", "index_prefix"]]
        ) | set{samples}

        sourcePrefix(
            BwaMemIndex,
            samples,
            "index_dir",
            "index_prefix",
            ["genomeReference", "index_prefix"],
            ["index_dir", "index_prefix", ".0123", ".ann", ".amb", ".pac", ".bwt", ".sa"],
            {["index_prefix":it.assembly]},
            "index_prefix",
            {it.datatype in ["fq", "fastq"] && it.aligner == "bwa"},
            ["keep":["index_dir", "index_prefix"]]
        ) | set{samples}

    samples = emptyOnLastStep("index", samples)

    emit:
        samples

}