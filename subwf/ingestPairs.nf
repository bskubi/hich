include {JoinProcessResults} from './joinProcessResults.nf'
include {QCReads} from './qcHicReads.nf'

process PairtoolsFlipSort {
    publishDir params.general.publish.flip_sort ? params.general.publish.flip_sort : "results",
               saveAs: {params.general.publish.flip_sort ? it : null}
    container "bskubi/pairtools:1.0.4"

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
        samples
            | filter{it.datatype == "pairs" && it.get("pairs") && file(it.get("pairs")).exists()}
            | map{
                it.pairs = file(it.pairs);
                it.id = it.sample_id;
                it
            }
            | set{ingest}
        
        samples = JoinProcessResults(
            PairtoolsFlipSort,
            [ingest, samples],
            ["sample_id", "pairs", "chromsizes"],
            ["sample_id", "pairs"],
            ["sample_id"],
            null,
            "pairs")
        
        if (params.general.get("last_step") == "IngestPairs") {
            channel.empty() | set{samples}
        }

        if ("IngestPairs" in params.general.get("qc_after")) {
            samples = QCReads(samples, "IngestPairs")
        }

    emit:
        samples
}