include {QCReads} from './qcHicReads.nf'
include {transpack} from './extraops.nf'

process PairtoolsParse2 {
    publishDir params.general.publish.parse ? params.general.publish.parse : "results",
               saveAs: {params.general.publish.parse ? it : null}
    container "bskubi/hich:latest"

    input:
    tuple val(sample_id), path(sambam), path(chromsizes), val(assembly), val(parse_params)

    output:
    tuple val(sample_id), path("${sample_id}.pairs.gz")

    shell:
    cmd = ["pairtools parse2",
           "--assembly ${assembly}",
           "--chroms-path ${chromsizes}"] +
           parse_params +
          ["${sambam} | pairtools sort --output ${sample_id}.pairs.gz"]
    cmd.removeAll([null])

    cmd.join(" ")

    stub:
    "touch ${sample_id}.pairs.gz"
}

workflow Parse {
    take:
    samples

    main:

    // This should give an error if the file does not exist
    samples
        | filter{it.get("sambam") && file(it.sambam).exists()}
        | map{it.sambam = file(it.sambam); it}
        | set {sambam}



    samples = transpack(
        PairtoolsParse2,
        [sambam, samples],
        ["sample_id", "sambam", "chromsizes", "assembly", "parse_params"],
        ["sample_id", "pairs"],
        ["latest":"pairs"],
        "sample_id"
        )
    
    samples | map{it.id = it.sample_id; it} | set{samples}

    // It might be good to simplify these workflow control steps since they
    // are repeated frequently.
    if ("Parse" in params.general.get("qc_after")) {
        QCReads(samples, "Parse")
    }

    if (params.general.get("last_step") == "parse") {
        channel.empty() | set{samples}
    }

    emit:
        samples
}
