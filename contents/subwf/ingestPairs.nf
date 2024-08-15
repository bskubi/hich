include {QCReads} from './qcHicReads.nf'
include {transpack; emptyOnLastStep} from './extraops.nf'

process PairtoolsFlipSort {
    publishDir params.general.publish.flip_sort ? params.general.publish.flip_sort : "results",
               saveAs: {params.general.publish.flip_sort ? it : null},
               mode: params.general.publish.mode
    conda "bioconda::pairtools"
    container "bskubi/hich:latest"

    input:
    tuple val(id), path(pairs), path(chromsizes)

    output:
    tuple val(id), path("${id}.pairs.gz")

    shell:
    cmd = ["pairtools flip --chroms-path ${chromsizes} ${pairs}",
           "| pairtools sort --output ${id}.pairs.gz"].join(" ")
    cmd
}

workflow IngestPairs {
    take:
        samples

    main:
    // Maybe there is an easier and more general way we can do ingestion?
    samples
        | filter{it.datatype == "pairs" && it.get("pairs") && file(it.get("pairs")).exists()}
        | map{
            it.pairs = file(it.pairs);
            it
        }
        | set{ingest}

    samples = transpack(
        PairtoolsFlipSort,
        [ingest, samples],
        ["id", "pairs", "chromsizes"],
        ["id", "pairs"],
        ["latest":"pairs"],
        "id"
    )

    if ("IngestPairs" in params.general.get("qc_after")) {
        samples = QCReads(samples, "IngestPairs")
    }

    samples = emptyOnLastStep("ingestPairs") ?: samples

    emit:
        samples
}