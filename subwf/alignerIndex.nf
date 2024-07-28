include {sqljoin} from './extraops.nf'

process BwaMem2Index {
    publishDir "resources/index/bwa-mem2/", mode: 'move'
    container "bskubi/bwa-mem2"

    input:
    tuple path(reference), val(prefix)

    output:
    tuple val(prefix), path("${prefix}.0123"), path("${prefix}.amb"),
          path("${prefix}.ann"), path("${prefix}.bwt.2bit.64"),
          path("${prefix}.pac")

    shell:
    "touch ${prefix}.0123 ${prefix}.amb ${prefix}.ann ${prefix}.bwt.2bit.64 ${prefix}.pac"
    //"bwa-mem2 index -p ${prefix} ${reference}"

    stub:
    "touch ${prefix}.0123 ${prefix}.amb ${prefix}.ann ${prefix}.bwt.2bit.64 ${prefix}.pac"
}

workflow MakeMissingIndex {

    take:
        samples
    
    main:

        def fileExists = { dir, file -> new File(dir, file).exists() }
        
        samples
            | filter{it.get("index_dir").getClass() != file("helloWorld").getClass()
                     && (it.get("index_dir", "").trim().length() == 0
                     || it.get("index_prefix").trim().length() == 0
                     || !fileExists(it.index_dir, "${it.index_prefix}.0123"))}
            | set{no_index}
        
        no_index
            | map{it.subMap("reference", "assembly")}
            | unique
            | map{it.index_prefix = it.index_prefix?.trim() ? it.index_prefix : it.assembly; it}
            | map{tuple(it.reference, it.index_prefix)}
            | BwaMem2Index
            | map{["index_prefix": it[0], "index_dir": "resources/index/bwa-mem2/"]}
            | set {new_index}

        sqljoin(no_index, new_index, [by: "assembly", suffix: ""])
            | set {new_index}

        sqljoin(samples, new_index, [by: "sample_id", suffix: ""])
            | map{it.index_dir = file(it.index_dir); it}
            | set {samples}

    emit:
        samples

}
