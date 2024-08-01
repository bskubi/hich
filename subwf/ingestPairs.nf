include {QCReads} from './qcHicReads.nf'
include {transpack} from './extraops.nf'

process PairtoolsFlipSort {
    publishDir params.general.publish.flip_sort ? params.general.publish.flip_sort : "results",
               saveAs: {params.general.publish.flip_sort ? it : null},
               mode: params.general.publish.mode
    conda "pairtools"
    container "bskubi/hich:latest"

    input:
    tuple val(sample_id), path(pairs), path(chromsizes)

    output:
    tuple val(sample_id), path("${sample_id}.pairs.gz")

    shell:
    cmd = ["pairtools flip --chroms-path ${chromsizes} ${pairs}",
           "| pairtools sort --output ${sample_id}.pairs.gz"]
    cmd.removeAll([null])
    cmd.join(" ")
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
                it.id = it.sample_id;
                it
            }
            | set{ingest}

        samples = transpack(
            PairtoolsFlipSort,
            [ingest, samples],
            ["sample_id", "pairs", "chromsizes"],
            ["sample_id", "pairs"],
            ["latest":"pairs"],
            "sample_id"
        )

        if ("IngestPairs" in params.general.get("qc_after")) {
            samples = QCReads(samples, "IngestPairs")
        }

        if (params.general.get("last_step") == "ingestpairs") {
            channel.empty() | set{samples}
        }

    emit:
        samples
}