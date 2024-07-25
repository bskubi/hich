include {JoinProcessResults} from './joinProcessResults.nf'
include {QCReads} from './qcHicReads.nf'

process PairtoolsParse2 {
    publishDir params.general.publish.parse ? params.general.publish.parse : "results",
               saveAs: {params.general.publish.parse ? it : null}
    container "bskubi/pairtools:1.1.0"
    //container "bskubi/pairtools:1.0.4"

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
}

workflow Parse {
    take:
    samples

    main:

        samples
        | filter{it.get("sambam") && file(it.sambam).exists()}
        | map{it.sambam = file(it.sambam); it}
        | set {sambam}
    
    samples = JoinProcessResults(
        PairtoolsParse2,
        [sambam, samples],
        ["sample_id", "sambam", "chromsizes", "assembly", "parse_params"],
        ["sample_id", "pairs"],
        ["sample_id"],
        null,
        "pairs")
    
    samples | map{it.id = it.sample_id; it} | set{samples}

    if ("Parse" in params.general.get("qc_after")) {
        QCReads(samples, "Parse")
    }

    if (params.general.get("last_step") == "parse") {
        channel.empty() | set{samples}
    }

    emit:
        samples
}
