include {transpack; emptyOnLastStep} from './extraops.nf'


process CoolerZoomify {
    publishDir params.general.publish.mcool ? params.general.publish.mcool : "results",
               saveAs: {params.general.publish.mcool ? it : null},
               mode: params.general.publish.mode
    conda "bioconda::cooler"
    container "bskubi/hich:latest"
    maxForks 2

    input:
    tuple val(id), path(infile), path(chromsizes), val(pairsFormat), val(assembly), val(matrix), val(coolerCloadParams), val(coolerZoomifyParams)

    output:
    tuple val(id), path("${id}.mcool")

    shell:
    min_bin = matrix.resolutions.min()
    bins = matrix.resolutions.join(",")

    cmd = ["cooler cload pairs"] + coolerCloadParams +
           ["--assembly ${assembly}",
           "--chrom1 ${pairsFormat.chrom1}",
           "--pos1 ${pairsFormat.pos1}",
           "--chrom2 ${pairsFormat.chrom2}",
           "--pos2 ${pairsFormat.pos2}",
           "${chromsizes}:${min_bin}",
           "${infile} ${id}.cool",
           "&& cooler zoomify"] + coolerZoomifyParams +
           ["--resolutions '${bins}'",
           "--out ${id}.mcool",
           "${id}.cool"]
    cmd.removeAll([null])
    cmd = cmd.join(" ")
    cmd

    stub:
    "touch ${id}.mcool"
}

workflow McoolMatrix {
    take:
        samples
    
    main:
        samples | filter{it.matrix.makeMcoolFileFormat} | set{mcool}

        transpack (
            CoolerZoomify,
            [mcool, samples],
            ["id", "latest", "chromsizes", "pairsFormat", "assembly", "matrix", "coolerCloadParams", "coolerZoomifyParams"],
            ["id", "mcool"],
            ["latestMatrix":"mcool"],
            "id",
            ["nullOk":"coolerCloadParams"]
        ) | set{samples}
    
    samples = emptyOnLastStep("McoolMatrix", samples)

    emit:
        samples
}