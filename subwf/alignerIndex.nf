include {sqljoin; sourcePrefix} from './extraops.nf'

process BwaMem2Index {
    container "bskubi/bwa-mem2"

    input:
    tuple path(reference), val(prefix)

    output:
    tuple path("bwa-mem2/index"), val(prefix), path("bwa-mem2/index/${prefix}.0123"), path("bwa-mem2/index/${prefix}.amb"),
          path("bwa-mem2/index/${prefix}.ann"), path("bwa-mem2/index/${prefix}.bwt.2bit.64"),
          path("bwa-mem2/index/${prefix}.pac")

    shell:
    "mkdir -p bwa-mem2/index && cd bwa-mem2/index && touch ${prefix}.0123 ${prefix}.amb ${prefix}.ann ${prefix}.bwt.2bit.64 ${prefix}.pac"
    //"bwa-mem2 index -p ${prefix} ${reference}"

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


// workflow MakeMissingIndex {

//     take:
//         samples
    
//     main:
//         /*
//             I reworked BwaMem2Index to create files in bwa-mem2/index/prefix.ext and to
//             return a path of the directory of bwa-mem2/index as well as a val of the prefix.

//             There's clear overlap with source, but a lot of changes have to be made.

//             Actually, the filter_unique step from source may not matter at all because
//             it will just use cached results. OK there is something weird going on with pairtools stats.
//             It seems to work fine except goes to hundreds of cached calls to pairtoolsStats if I remove
//             the unique section.

//             needs = {it.datatype == "fastq"}
//             namer = {[index_dir:"bwa-mem2/index", index_prefix:it.assembly]}

//             We need to produce two attributes:
//             index_dir and index_prefix.

//             Rework namer to return a hashmap update instead of values for the key
//             Pass a caster closure that casts some hashmap values to files.
            
//             Set as missing if the type is XPath or is a file that does not exist







//         */


//         def fileExists = { dir, file -> new File(dir, file).exists() }
        
//         samples
//             | filter{it.get("index_dir").getClass() != file("helloWorld").getClass()
//                      && (it.get("index_dir", "").trim().length() == 0
//                      || it.get("index_prefix").trim().length() == 0
//                      || !fileExists(it.index_dir, "${it.index_prefix}.0123"))}
//             | set{no_index}
        
//         no_index
//             | map{it.subMap("reference", "assembly")}
//             | unique
//             | map{it.index_prefix = it.index_prefix?.trim() ? it.index_prefix : it.assembly; it}
//             | map{tuple(it.reference, it.index_prefix)}
//             | BwaMem2Index
//             | map{["index_prefix": it[0], "index_dir": "resources/index/bwa-mem2/"]}
//             | set {new_index}

//         sqljoin(no_index, new_index, [by: "assembly", suffix: ""])
//             | set {new_index}

//         sqljoin(samples, new_index, [by: "sample_id", suffix: ""])
//             | map{it.index_dir = file(it.index_dir); it}
//             | set {samples}

//     emit:
//         samples

// }
