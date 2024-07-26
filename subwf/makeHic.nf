include {JoinProcessResults} from './joinProcessResults.nf'
include {transpack; hashmapdiff} from './extraops.nf'

process JuicerToolsPre {
    /*
        Juicer Tools Pre documentation: https://github.com/aidenlab/juicer/wiki/Pre

        We use version 1.22.01 as 2.0+ versions are in development and certain
        features available in version 1 are unavailable in 2.
    */
    publishDir params.general.publish.hic ? params.general.publish.hic : "results",
               saveAs: {params.general.publish.hic ? it : null}
    container "bskubi/juicer_tools:1.22.01"
    maxForks 2

    input:
    tuple val(id), path(infile), path(chromsizes), val(pairs_format), val(matrix)

    output:
    tuple val(id), path("${id}.hic")

    shell:
    outfile = "${id}.hic"
    min_bin = matrix.resolutions.min()
    bins = matrix.resolutions ? "-r ${matrix.resolutions.join(',')}" : ""

    cmd = ["java -Xmx20g -jar /app/juicer_tools.jar pre",
            bins,
           "${infile} ${outfile} ${chromsizes}"]
    cmd.removeAll([null])
    cmd.join(" ")

    stub:
    "touch ${id}.hic"
}

workflow MakeHic {
    take:
        samples
    
    main:
        samples | filter{it.matrix.make_hic_file_format} | set{hic}

        transpack (
            JuicerToolsPre,
            [hic, samples],
            ["id", "latest", "chromsizes", "pairs_format", "matrix"],
            ["id", "hic"],
            ["latest_matrix":"hic"]
        ) | set{samples}

    emit:
        samples
}