include {QCReads} from './qcHicReads.nf'
include {transpack; emptyOnLastStep} from './extraops.nf'

process PairtoolsDedup {
    publishDir params.general.publish.dedup ? params.general.publish.dedup : "results",
               saveAs: {params.general.publish.dedup ? it : null},
               mode: params.general.publish.mode
    conda "bioconda::pairtools"
    container "bskubi/hich:latest"

    input:
    tuple val(id), path(infile), val(pairtoolsDedupParams)

    output:
    tuple val(id), path("${id}_dedup.pairs.gz")

    shell:
    pairtoolsDedupParams = pairtoolsDedupParams ? pairtoolsDedupParams.collect
        {
            item ->
            return [
                "--output-stats": "--output-stats ${id}_dedup.stats.txt",
                "--output-dups": "--output-dups ${id}_dedup.dups.pairs.gz",
                "--output-unmapped": "--output-unmapped ${id}_dedup.unmapped.pairs.gz",
                "--output-bytile-stats": "--output-bytile-stats ${id}_dedup.bytile_stats.pairs.gz"
            ].get(item, item)
        }.join(" ") : ""
    
    cmd = "pairtools dedup --output ${id}_dedup.pairs.gz ${pairtoolsDedupParams} ${infile}"
    cmd

    stub:
    "touch ${id}_dedup.pairs.gz"
}

workflow Deduplicate {
    take:
        samples
    
    main:
        
    samples | filter{it.deduplicate} | set {deduplicate}

    samples = transpack(
        PairtoolsDedup,
        [deduplicate, samples],
        ["id", "latest", "pairtoolsDedupParams"],
        ["id", "dedup_pairs"],
        ["latest":"dedup_pairs"],
        "id",
        ["nullOk":"pairtoolsDedupParams"]
    )
    
    samples = emptyOnLastStep("Deduplicate", samples)

    emit:
        samples
}