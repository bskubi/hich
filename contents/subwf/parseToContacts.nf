include {QCReads} from './qcHicReads.nf'
include {transpack} from './extraops.nf'

process PairtoolsParse2 {
    publishDir params.general.publish.parse ? params.general.publish.parse : "results",
               saveAs: {params.general.publish.parse ? it : null},
               mode: params.general.publish.mode
    conda "bioconda::pairtools"
    container "bskubi/hich:latest"

    input:
    tuple val(id), path(sambam), path(chromsizes), val(assembly), val(parse_params)

    output:
    tuple val(id), path("${id}.pairs.gz")

    shell:
    
    cmd = ["samtools sort -n ${sambam}",
           "| pairtools parse2",
           "--assembly ${assembly}",
           "--chroms-path ${chromsizes}"] +
           parse_params +
          ["| pairtools sort --output ${id}.pairs.gz"]
    cmd.removeAll([null])

    cmd.join(" ")

    stub:
    "touch ${id}.pairs.gz"
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
        ["id", "sambam", "chromsizes", "assembly", "parse_params"],
        ["id", "pairs"],
        ["latest":"pairs"],
        "id"
        )

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
