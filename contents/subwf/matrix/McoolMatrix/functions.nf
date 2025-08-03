def buildCmd(id, matrixPlanName, pairs, chromsizes, assembly, matrix, coolerCloadParams, coolerZoomifyParams, cpus) {
    coolerCloadParams = coolerCloadParams ?: []
    coolerZoomifyParams = coolerZoomifyParams ?: []
    def resolutions = matrix.resolutions

    // Extract resolutions from --resolutions or -r parameter if specified
    // Otherwise use resolutions
    def bins = []
    def resolutionsFromParams = coolerZoomifyParams.find{it.startsWith("--resolutions") || it.startsWith("-r")}
    if (resolutionsFromParams) {
        try {
            bins = resolutionsFromParams.split()[1].split(",").collect{it.toInteger()}
        } catch(NumberFormatException e) {
            print("In CoolerZoomify on id ${id}, resolutions was ${resolutionsFromParams}, which needs to be a comma-separated list of integer resolutions")
        }        
    } else {
        bins = resolutions
        coolerZoomifyParams += ["--resolutions ${bins.join(",")}"]
    }
    

    if (!coolerZoomifyParams.any{it.startsWith("--nproc") || it.startsWith("-n") || it.startsWith("-p")}) {
        coolerZoomifyParams += ["--nproc ${cpus}"]
    }

    // Default to using assembly for --assembly argument
    if (!coolerCloadParams.any{it.startsWith("--assembly")} && assembly) {
        coolerCloadParams += ["--assembly ${assembly}"] 
    }

    assert bins, "In CoolerZoomify on id ${id}, matrix.resolutions is unspecified and no --resolutions is given."
    assert bins.every{it instanceof Integer}, "In CoolerZoomify on id ${id}, resolutions ${resolutions}, matrix.resolutions ${matrix.resolutions}, bins was parsed as ${bins} which contains a non-integer"
    def minBin = bins.min()
    def cmd = (["cooler cload pairs"]
           + coolerCloadParams
           + ["'${chromsizes}:${minBin}'", "'${pairs}' '${id}.cool'", "&& cooler zoomify"]
           + coolerZoomifyParams
           + ["--out '${id}.mcool'", "'${id}.cool'"])
    cmd = cmd.findAll{it}
    cmd = cmd.join(" ")
    output = "${id}.mcool"

    logMap = [task: "MCOOL_MATRIX", output: output, input: [id: id, pairs: pairs, chromsizes: chromsizes, assembly: assembly, matrix: matrix, coolerCloadParams: coolerCloadParams, coolerZoomifyParams: coolerZoomifyParams]]
    

    return [cmd, logMap, output]
}