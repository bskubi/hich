include {sqljoin; sourcePrefix} from './extraops.nf'

process BwaMem2Index {
    container "bskubi/bwa-mem2"
    publishDir params.general.publish.index ? params.general.publish.index : "results",
               saveAs: {params.general.publish.index ? it : null}

    input:
    tuple path(reference), val(prefix)

    output:
    tuple path("bwa-mem2/index"), val(prefix), path("bwa-mem2/index/${prefix}.0123"), path("bwa-mem2/index/${prefix}.amb"),
          path("bwa-mem2/index/${prefix}.ann"), path("bwa-mem2/index/${prefix}.bwt.2bit.64"),
          path("bwa-mem2/index/${prefix}.pac")

    shell:
    //"mkdir -p bwa-mem2/index && cd bwa-mem2/index && touch ${prefix}.0123 ${prefix}.amb ${prefix}.ann ${prefix}.bwt.2bit.64 ${prefix}.pac"
    "bwa-mem2 index -p ${prefix} ${reference}"

    stub:
    "mkdir -p bwa-mem2/index && cd bwa-mem2/index && touch ${prefix}.0123 ${prefix}.amb ${prefix}.ann ${prefix}.bwt.2bit.64 ${prefix}.pac"
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