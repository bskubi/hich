include {QCReads} from './qcHicReads.nf'
include {transpack; emptyOnLastStep; isSingleCell} from './extraops.nf'

process PairtoolsDedup {
    publishDir params.general.publish.dedup ? params.general.publish.dedup : "results",
               saveAs: {params.general.publish.dedup ? it : null},
               mode: params.general.publish.mode
    conda "bioconda::pairtools"
    container "bskubi/hich:latest"

    input:
    tuple val(id), path(infile), val(pairtoolsDedupParams), val(isSingleCell)

    output:
    tuple val(id), path("${id}_dedup.pairs.gz")

    shell:
    if (isSingleCell) {
        /*
            The --extra-col-pair cellID cellID option enforces that reads must have
            the same cellID to be called as duplicates.

            --backend cython makes single-cell deduplication non-transitive and
            performant. The default backend is far too slow.
        */
        pairtoolsDedupParams += ["--extra-col-pair cellID cellID --backend cython --chunksize 1000000"]
    }

    /*
        There are several output files for stats, dups, unmapped, and bytile_stats
        that can optionally be included in the dedup output.

        This section builds the filenames for any of these that the user requests
        based on the id, just like with the main output filename, leaving
        other options unaltered.
    */
    pairtoolsDedupParams = pairtoolsDedupParams ? pairtoolsDedupParams.collect
        {
            item ->

            // Item is a single param. The following is a HashMap of keys
            // linking the flag to the flag + id-based filename. If 'item' is one of the
            // keys, then it returns the flag + id-based filename. Otherwise,
            // it returns the item itself.
            return [
                "--output-stats": "--output-stats ${id}_dedup.stats.txt",
                "--output-dups": "--output-dups ${id}_dedup.dups.pairs.gz",
                "--output-unmapped": "--output-unmapped ${id}_dedup.unmapped.pairs.gz",
                "--output-bytile-stats": "--output-bytile-stats ${id}_dedup.bytile_stats.pairs.gz"
            ].get(item, item)
        }.join(" ") : ""
    
    // Build the full command and run
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
        ["id", "latest", "pairtoolsDedupParams", "isSingleCell"],
        ["id", "dedup_pairs"],
        ["latest":"dedup_pairs"],
        "id",
        ["nullOk":["pairtoolsDedupParams", "isSingleCell"]]
    )
    
    samples = emptyOnLastStep("Deduplicate", samples)

    emit:
        samples
}