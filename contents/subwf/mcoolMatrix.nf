include {emptyOnLastStep; pack; skip} from './extraops.nf'


process CoolerZoomify {
    publishDir params.general.publish.mcool ? params.general.publish.mcool : "results",
               saveAs: {params.general.publish.mcool ? it : null},
               mode: params.general.publish.mode
    conda "bioconda::cooler"
    container params.general.hichContainer
    label 'createMatrix'

    input:
    tuple val(id), path(infile), path(chromsizes), val(pairsFormat), val(assembly), val(matrix), val(coolerCloadParams), val(coolerZoomifyParams)

    output:
    tuple val(id), path("${id}.mcool")

    shell:
    coolerCloadParams = coolerCloadParams ?: []
    coolerZoomifyParams = coolerZoomifyParams ?: []

    // Extract resolutions from --resolutions or -r parameter if specified
    // Otherwise use matrix.resolutions
    bins = []
    resolutions = coolerZoomifyParams.find{it.startsWith("--resolutions") || it.startsWith("-r")}
    if (resolutions) {
        try {
            bins = resolutions.split()[1].split(",").collect{it.toInteger()}
        } catch(NumberFormatException e) {
            print("In CoolerZoomify on id ${id}, resolutions was ${resolutions}, which needs to be a comma-separated list of integer resolutions")
        }        
    } else {
        bins = matrix.resolutions
        coolerZoomifyParams += ["--resolutions ${bins.join(",")}"]
    }
    

    // Default to using task.cpus as number of processes
    if (!coolerZoomifyParams.any{it.startsWith("--nproc") || it.startsWith("-n") || it.startsWith("-p")}) {
        coolerZoomifyParams += ["--nproc ${task.cpus}"]
    }

    // Default to using assembly for --assembly argument
    if (!coolerCloadParams.any{it.startsWith("--assembly")} && assembly) {
        coolerCloadParams += ["--assembly ${assembly}"] 
    }

    assert bins, "In CoolerZoomify on id ${id}, matrix.resolutions is unspecified and no --resolutions is given."
    assert bins.every{it instanceof Integer}, "In CoolerZoomify on id ${id}, resolutions ${resolutions}, matrix.resolutions ${matrix.resolutions}, bins was parsed as ${bins} which contains a non-integer"
    
    min_bin = bins.min()

    cmd = (["cooler cload pairs"]
           + coolerCloadParams
           + ["'${chromsizes}:${min_bin}'", "'${infile}' '${id}.cool'", "&& cooler zoomify"]
           + coolerZoomifyParams
           + ["--out '${id}.mcool'", "'${id}.cool'"])
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
    samples
    | filter{!skip("mcoolMatrix") && it.matrix.makeMcoolFileFormat && (it.pairs || it.latestPairs) && !it.mcool}
    | map{tuple(it.id, it.latest, it.chromsizes, it.pairsFormat, it.assembly, it.matrix, it.coolerCloadParams, it.coolerZoomifyParams)}
    | CoolerZoomify
    | map{id, mcool -> [id: id, mcool: mcool, latestMatrix: mcool]}
    | set{result}
    pack(samples, result) | set{samples}
    
    samples = emptyOnLastStep("mcoolMatrix", samples)

    emit:
    samples
}