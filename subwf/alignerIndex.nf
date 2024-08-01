include {sqljoin; sourcePrefix} from './extraops.nf'

process BwaMem2Index {
    container "bskubi/hich:latest"
    publishDir params.general.publish.bwa_mem2_index ? params.general.publish.bwa_mem2_index : "results",
               saveAs: {params.general.publish.bwa_mem2_index ? it : null},
               mode: params.general.publish.mode

    input:
    tuple path(reference), val(prefix)

    output:
    tuple path("bwa-mem2/index"), val(prefix), path("bwa-mem2/index/${prefix}.0123"), path("bwa-mem2/index/${prefix}.amb"),
          path("bwa-mem2/index/${prefix}.ann"), path("bwa-mem2/index/${prefix}.bwt.2bit.64"),
          path("bwa-mem2/index/${prefix}.pac")

    shell:
    //"mkdir -p bwa-mem2/index && cd bwa-mem2/index && touch ${prefix}.0123 ${prefix}.amb ${prefix}.ann ${prefix}.bwt.2bit.64 ${prefix}.pac"
    "bwa-mem2 index -p ${prefix} ${reference} && mkdir -p bwa-mem2/index && mv ${prefix}.0123 ${prefix}.amb ${prefix}.ann ${prefix}.bwt.2bit.64 ${prefix}.pac bwa-mem2/index"

    stub:
    "mkdir -p bwa-mem2/index && cd bwa-mem2/index && touch ${prefix}.0123 ${prefix}.amb ${prefix}.ann ${prefix}.bwt.2bit.64 ${prefix}.pac"
}

process BwaMemIndex {
    container "bskubi/hich:latest"
    publishDir params.general.publish.bwa_mem_index ? params.general.publish.bwa_mem2_index : "results",
               saveAs: {params.general.publish.bwa_mem_index ? it : null},
               mode: params.general.publish.mode

    input:
    tuple path(reference), val(prefix)

    output:
    tuple path("bwa/index"), val(prefix), path("bwa/index/${prefix}.ann"), path("bwa/index/${prefix}.amb"),
          path("bwa/index/${prefix}.pac"), path("bwa/index/${prefix}.bwt"), path("bwa/index/${prefix}.sa"),

    shell:
    //"mkdir -p bwa-mem2/index && cd bwa-mem2/index && touch ${prefix}.0123 ${prefix}.amb ${prefix}.ann ${prefix}.bwt.2bit.64 ${prefix}.pac"
    "bwa index -p ${prefix} ${reference} && mkdir -p bwa/index && mv ${prefix}.amb ${prefix}.ann ${prefix}.pac ${prefix}.bwt ${prefix}.sa bwa/index"

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
            ["reference", "index_prefix"],
            ["index_dir", "index_prefix", ".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac"],
            {["index_prefix":it.assembly]},
            "index_prefix",
            {it.datatype in ["fq", "fastq"]},
            ["keep":["index_dir", "index_prefix"]]
        ) | set{samples}

    emit:
        samples

}