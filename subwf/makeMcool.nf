include {transpack} from './extraops.nf'


process CoolerZoomify {
    publishDir params.general.publish.mcool ? params.general.publish.mcool : "results",
               saveAs: {params.general.publish.mcool ? it : null}
    conda "cooler"
    maxForks 2

    input:
    tuple val(id), path(infile), path(chromsizes), val(pairs_format), val(assembly), val(matrix), val(make_cool), val(make_mcool)

    output:
    tuple val(id), path("${id}.mcool")

    shell:
    min_bin = matrix.resolutions.min()
    bins = matrix.resolutions.join(",")

    cmd = ["cooler cload pairs"] + make_cool +
           ["--assembly ${assembly}",
           "--chrom1 ${pairs_format.chrom1}",
           "--pos1 ${pairs_format.pos1}",
           "--chrom2 ${pairs_format.chrom2}",
           "--pos2 ${pairs_format.pos2}",
           "${chromsizes}:${min_bin}",
           "${infile} ${id}.cool",
           "&& cooler zoomify"] + make_mcool +
           ["--resolutions '${bins}'",
           "--out ${id}.mcool",
           "${id}.cool"]
    cmd.removeAll([null])
    cmd = cmd.join(" ")
    cmd

    stub:
    "touch ${id}.mcool"
}

workflow MakeMcool {
    take:
        samples
    
    main:
        samples
            | filter{it.matrix.make_mcool_file_format}
            | set{mcool}

        transpack (
            CoolerZoomify,
            [mcool, samples],
            ["id", "latest", "chromsizes", "pairs_format", "assembly", "matrix", "make_cool", "make_mcool"],
            ["id", "mcool"],
            ["latest_matrix":"mcool"],
            "id",
            ["nullOk":"make_cool"]
        ) | set{samples}

    emit:
        samples
}