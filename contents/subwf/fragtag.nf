include {QCReads} from './qcHicReads.nf'
include {transpack; emptyOnLastStep} from './extraops.nf'

process FragtagProc {
    publishDir params.general.publish.fragtag ? params.general.publish.fragtag : "results",
               saveAs: {params.general.publish.fragtag ? it : null},
               mode: params.general.publish.mode
    //container "bskubi/pairtools:1.1.0"
    container "bskubi/hich:latest"

    input:
    tuple val(id), path(pairs), path(fragmentIndex)

    output:
    tuple val(id), path("${id}_fragtag.pairs.gz")

    shell:
    // ["pairtools restrict",
    //  "--frags ${fragfile}",
    //  "--output ${tagged_pairs}",
    //  "${pairs}"].join(" ")
    cmd = "hich fragtag ${fragmentIndex} ${id}_fragtag.pairs.gz ${pairs}"
    cmd

    stub:
    "touch ${id}_fragtag.pairs.gz"
}

def hasFragmentIndexName = {it.get("fragmentIndex").toString().trim().length() > 0}

def fragmentIndexExists = {hasFragmentIndexName(it) && file(it.fragmentIndex).exists()}

workflow Fragtag {
    take:
        samples

    main:
    samples | filter{fragmentIndexExists(it)} | set{fragtag}

    samples = transpack(
        FragtagProc,
        [fragtag, samples],
        ["id", "pairs", "fragmentIndex"],
        ["id", "frag_pairs"],
        ["latest":"frag_pairs"],
        "id"
    )

    if ("OptionalFragtag" in params.general.get("qc_after")) {
        samples = QCReads(samples, "OptionalFragtag")
    }

    samples = emptyOnLastStep("fragtag", samples)

    emit:
        samples
}