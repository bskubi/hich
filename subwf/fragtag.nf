include {QCReads} from './qcHicReads.nf'
include {transpack} from './extraops.nf'

process Fragtag {
    publishDir params.general.publish.fragtag ? params.general.publish.fragtag : "results",
               saveAs: {params.general.publish.fragtag ? it : null},
               mode: params.general.publish.mode
    //container "bskubi/pairtools:1.1.0"
    container "bskubi/hich:latest"

    input:
    tuple val(sample_id), path(pairs), path(fragfile), val(tagged_pairs)

    output:
    tuple val(sample_id), path(tagged_pairs)

    shell:
    // ["pairtools restrict",
    //  "--frags ${fragfile}",
    //  "--output ${tagged_pairs}",
    //  "${pairs}"].join(" ")
    
    cmd = "hich fragtag ${fragfile} ${tagged_pairs} ${pairs}"
    cmd

    stub:
    "touch ${tagged_pairs}"
}

workflow OptionalFragtag {
    take:
        samples

    main:
        // I think this is cruft and can be removed
        def hasFragfileName = {
            it.get("fragfile").toString().trim().length() > 0
        }
        
        def fragfileExists = {
            hasFragfileName(it) && file(it.fragfile).exists()
        }
        // End cruft

        samples
            | filter{fragfileExists(it)}
            | map{it.frag_pairs = "${it.sample_id}_fragtag.pairs.gz"; it}
            | set{fragtag}

        samples = transpack(
            Fragtag,
            [fragtag, samples],
            ["sample_id", "pairs", "fragfile", "frag_pairs"],
            ["sample_id", "frag_pairs"],
            ["latest":"frag_pairs"],
            "sample_id"
        )

        if ("OptionalFragtag" in params.general.get("qc_after")) {
            samples = QCReads(samples, "OptionalFragtag")
        }

    emit:
        samples
}