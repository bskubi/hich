include {QCReads} from './qcHicReads.nf'
include {transpack; emptyOnLastStep} from './extraops.nf'

process PairtoolsParse2 {
    publishDir params.general.publish.parse ? params.general.publish.parse : "results",
               saveAs: {params.general.publish.parse ? it : null},
               mode: params.general.publish.mode
    conda "bioconda::pairtools bioconda::samtools"
    container "bskubi/hich:latest"

    input:
    tuple val(id), path(sambam), path(chromsizes), val(assembly), val(parseParams), val(reshapeParams)

    output:
    tuple val(id), path("${id}.pairs.gz")

    shell:
    // Previous implementation
    // cmd = ["samtools sort -n ${sambam}",    
    //        "| pairtools parse2",            
    //        "--assembly ${assembly}",
    //        "--chroms-path ${chromsizes}"] +
    //        parseParams +                               
    //       ["| pairtools sort --output ${id}.pairs.gz"] 
    // cmd.removeAll([null])
    // cmd.join(" ")

    /*
        We combine several steps into one here for conceptual simplicity and
        to save on space by avoiding unnecessary duplication of the .pairs file
        for each step.

        Samtools sort sorts the reads by name so that paired alignments are grouped
        together. Otherwise, corrupt and inaccurate parsing results. This should
        be replaced by sambamba sort, which is faster. It takes place here rather
        than after alignment because ingested .sam/.bam files also need to have
        this sorting applied.

        pairtools parse2 parses the reads and rescues reads with observed ligation
        junctions (i.e. mapping to locus 1 AND 2 on one side, and to locus 1 OR 2 on the other).

        If the user wants to drop the readID column for single cell data, this takes place during
        hich reshape rather than pairtools parse2 in order to ensure the readID information is available.

        hich reshape is used for single cell data to extract a cell barcode from a given cell
        using a regex or Python parse library pattern and move it to the cellID column.
        This permits single-cell deduplication, where reads are not considered duplicates unless
        they have the same cellID (i.e. originate from the same cell). It also preserves the cell
        ID information through the rest of the pipeline so that users can perform single cell
        analysis.

        Per the pairtools spec, pairtools sort sorts by "lexicographic order along chrom1 and chrom2,
        in the numeric order along pos1 and pos2 and in the lexicographic order along pair_type."
    */

    // The parseParams are typically a list of individual pairtools parse flags.
    // Join them separated by spaces to use in the parse2Cmd
    parseParams = parseParams.join(" ")
    reshapeParams = reshapeParams.join(" ")

    // Set up the individual commands in lists to make them easier to combine with pipes into a complete command
    // sambamba is both slower than samtools as of 2017, and also can't pipe to stdout, so we use samtools
    samSortCmd = ["samtools sort -n ${sambam}"]
    parse2Cmd = ["pairtools parse2 --assembly ${assembly} --chroms-path ${chromsizes} ${parseParams}"]
    reshapeCmd = reshapeParams ? ["hich reshape ${reshapeParams}"] : []
    pairsSortCmd = ["pairtools sort --output ${id}.pairs.gz"]

    // Combine the individual commands, then join with a pipe to form the full command
    cmdParts = samSortCmd + parse2Cmd + reshapeCmd + pairsSortCmd
    cmd = cmdParts.join(" | ")
    

    // Execute the full command
    cmd

    stub:
    "touch ${id}.pairs.gz"
}

workflow Parse {
    take:
    samples

    main:
    // This should give an error if the file does not exist
    samples
        | filter{it.datatype == "sambam"}
        | map{
            sample ->
            if (!sample.sambam || !file(sample.sambam).exists()) {
                error "In sample with id ${sample.id}, sambam file is specified but does not exist"
            }

            sample
        }
        | set {sambam}



    samples = transpack(
        PairtoolsParse2,
        [sambam, samples],
        ["id", "sambam", "chromsizes", "assembly", "parseParams", "reshapeParams"],
        ["id", "pairs"],
        ["latest":"pairs"],
        "id",
        ["nullOk":["reshapeParams", "parseParams"]]
        )

    // It might be good to simplify these workflow control steps since they
    // are repeated frequently.
    if ("Parse" in params.general.get("qcAfter")) {
        QCReads(samples, "Parse")
    }

    samples = emptyOnLastStep("Parse", samples)

    emit:
        samples
}
