include {buildCLIOpts} from '../../util/cli.nf'

def buildCmdError(id, matrix_opts, message) {
    message = "Error in cooler zoomify resolution specification for id ${id}, matrix_opts ${matrix_opts}: ${message}"
    error(message)
}

def buildCmd(id, pairs, chromsizes, assembly, matrix_opts, cpus) {
    def mcool_output = "${id}.mcool"
    def logMap = [
        task: "MCOOL_MATRIX", 
        output: [mcool: mcool_output], 
        input: [
            id: id, 
            pairs: pairs, 
            chromsizes: chromsizes, 
            assembly: assembly,
            matrix_opts: matrix_opts
        ]
    ]

    def cooler_zoomify_opts = matrix_opts?.cooler_zoomify_opts ?: [:]

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
        "-p2": 5,
        "--temp-dir": "."
    ]
    def cooler_cload_pairs_opts = matrix_opts?.cooler_cload_pairs ?: [:]
    def bestResolution = resolutions.min()
    def chromsizesBestResolution = "${chromsizes}:${bestResolution}"
    def remap = [
        "--chrom1": "-c1",
        "--pos1": "-p1",
        "--chrom2": "-c2",
        "--pos2": "-p2",
        "--zero-based": "-0",
        "--no-symmetric-upper": "-N"
    ]
    logMap.default_cooler_cload_pairs_opts = default_cooler_cload_pairs_opts
    logMap.cooler_cload_pairs_opts = cooler_cload_pairs_opts
    cooler_cload_pairs_opts = buildCLIOpts(default_cooler_cload_pairs_opts, cooler_cload_pairs_opts, remap, null)
    def coolOutput = "${id}.cool"
    def args = [chromsizesBestResolution, pairs, coolOutput].collect{"'${it}'"}.join(" ")
    def cooler_cload_pairs_cmd = "cooler cload pairs ${cooler_cload_pairs_opts} ${args}"

    def default_cooler_balance_opts = [
        "--nproc": cpus,
        "--max-iters": 2000,
        "--trans-only": true
    ]
    def cooler_balance_opts = matrix_opts?.cooler_balance_opts ?: [:]
    remap = [
        "--chunksize": "-c",
        "--force": "-f"
    ]
    cooler_balance_opts = buildCLIOpts(default_cooler_balance_opts, cooler_balance_opts, remap, null)

    resolutions = resolutions.join(",")
    def default_cooler_zoomify_opts = [
        "--nproc": cpus,
        "--balance": true,
        "--balance-args": cooler_balance_opts,
        "--resolutions": resolutions,
        "--out": mcool_output
    ]
    reshape = [
        "--nproc": "-n",
        "-p": "-n",
        "--chunksize": "-c",
        "--resolutions": "-r"
    ]
    cooler_zoomify_opts = buildCLIOpts(default_cooler_zoomify_opts, cooler_zoomify_opts, reshape, ["--balance-args": "\""])
    def cooler_zoomify_cmd = "cooler zoomify ${cooler_zoomify_opts} '${coolOutput}'"

    def cmd = "${cooler_cload_pairs_cmd} && ${cooler_zoomify_cmd}"
    logMap.cmd = cmd
    return [cmd, logMap, mcool_output]
}