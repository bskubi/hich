include {buildCLIOpts} from '../../util/cli.nf'

def buildCmdError(id, matrixPlanName, matrix_opts, message) {
    message = "Error in cooler zoomify resolution specification for id ${id}, matrixPlanName ${matrixPlanName}, matrix_opts: ${matrix_opts}: ${message}"
    error(message)
}

def buildCmd(id, matrixPlanName, pairs, chromsizes, assembly, matrix_opts, cpus) {
    def cooler_zoomify_opts = matrix_opts?.cooler_zoomify ?: [:]
    if (!matrix_opts?.resolutions && !cooler_zoomify_opts["-r"] && !cooler_zoomify_opts["--resolutions"]) {
        buildCmdError(id, matrixPlanName, matrix_opts, "No resolutions specified for matrix_opts: ${matrix_opts}.")
    }
    def resolutions = matrix_opts?.resolutions ?: [100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000, 500000, 1000000, 2000000, 5000000]
    if ("--resolutions" in cooler_zoomify_opts) {
        resolutions = cooler_zoomify_opts["--resolutions"]
    }
    if ("-r" in cooler_zoomify_opts) {
        if (resolutions && cooler_zoomify_opts["-r"] != resolutions) {
            buildCmdError(id, matrixPlanName, matrix_opts, "matrix_opts.cooler_zoomify included both -r and --resolutions, with mismatched values.")
        } else {
            resolutions = cooler_zoomify_opts["-r"]
        }
    }

    def default_cooler_cload_pairs_opts = [
        "--assembly": assembly,
        "-c1": 2,
        "-p1": 3,
        "-c2": 4,
        "-p2": 5
    ]
    def cooler_cload_pairs_opts = matrix_opts?.cooler_cload_pairs ?: [:]
    def bestResolution = resolutions.min()
    def chromsizesBestResolution = "${chromsizes}:${bestResolution}"
    cooler_cload_pairs_opts = buildCLIOpts(default_cooler_cload_pairs_opts, cooler_cload_pairs_opts, null, null)
    def coolOutput = "${id}.cool"
    def coolerCloadPairsCmd = "cooler cload pairs ${cooler_cload_pairs_opts} '${chromsizesBestResolution}' '${pairs}' '${coolOutput}'"

    def default_cooler_balance_opts = [
        "--nproc": cpus,
        "--max-iters": 2000,
        "--trans-only": true
    ]
    def cooler_balance_opts = matrix_opts?.cooler_balance ?: [:]
    cooler_balance_opts = buildCLIOpts(default_cooler_balance_opts, cooler_balance_opts, null, null)
    resolutions = resolutions.join(",")
    def mcoolOutput = "${id}.mcool"
    def default_cooler_zoomify_opts = [
        "--nproc": cpus,
        "--balance": true,
        "--balance-args": cooler_balance_opts,
        "--resolutions": resolutions,
        "--out": mcoolOutput
    ]
    cooler_zoomify_opts = cooler_zoomify_opts.findAll{argName, argVal -> !(argName in  ["-r", "--resolutions"])}
    cooler_zoomify_opts = buildCLIOpts(default_cooler_zoomify_opts, cooler_zoomify_opts, ["--resolutions": "-r"], ["--balance-args": "\""])
    def coolerZoomifyCmd = "cooler zoomify ${cooler_zoomify_opts} '${id}.cool'"
    def cmd = "${coolerCloadPairsCmd} && ${coolerZoomifyCmd}"

    logMap = [
        task: "MCOOL_MATRIX", 
        output: mcoolOutput, 
        input: [
            id: id, 
            pairs: pairs, 
            chromsizes: chromsizes, 
            assembly: assembly,
            matrix_opts: matrix_opts
        ]
    ]

    return [cmd, logMap, mcoolOutput]
}