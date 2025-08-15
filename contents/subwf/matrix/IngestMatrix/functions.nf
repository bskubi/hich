include {buildCLIOpts} from '../../util/cli.nf'

def buildCmdMcoolToHic(id, mcoolFile, ingest_matrix_opts, cpus) {
    def output = "${id}.hic"
    def logMap = [
        task: "MCOOL_TO_HIC", 
        output: [hic: output], 
        input: [
            id: id, 
            mcoolFile: mcoolFile,
            ingest_matrix_opts: ingest_matrix_opts,
            cpus: cpus
        ],
        cpus: cpus
    ]
    
    cpus = cpus >= 2 ? cpus : 2

    def default_hictk_convert_mcool_to_hic_opts = [
        "-t": cpus,
        "--tmpdir":"."
    ]
    def hictk_convert_mcool_to_hic_opts = ingest_matrix_opts?.hictk_convert_mcool_to_hic_opts ?: [:]
    def remap = [
        "--resolutions": "-r",
        "--genome": "-g",
        "--verbosity": "-v",
        "--threads": "-t",
        "--compression-lvl": "-l",
        "--force": "-f"
    ]
    def final_hictk_convert_mcool_to_hic_opts = buildCLIOpts(default_hictk_convert_mcool_to_hic_opts, hictk_convert_mcool_to_hic_opts, remap, null)
    logMap += [default_hictk_convert_mcool_to_hic_opts: default_hictk_convert_mcool_to_hic_opts, hictk_convert_mcool_to_hic_opts: hictk_convert_mcool_to_hic_opts]
    def args = [mcoolFile, output].collect{"'${it}'"}.join(" ")
    def cmd = "hictk convert ${final_hictk_convert_mcool_to_hic_opts} ${args}"
    logMap += [cmd: cmd]
    return [cmd, logMap, output]
}

def buildCmdHicToMcool(id, hicFile, ingest_matrix_opts) {
    // Bugs result when trying to use >2 cpus for hictk convert from .hic to .mcool so we don't set this parameter
    // unless explicitly specified.
    def output = "${id}.mcool"
    def logMap = [
        task: "HIC_TO_MCOOL", 
        output: [mcool: output], 
        input: [
            id: id, 
            hicFile: hicFile,
            ingest_matrix_opts: ingest_matrix_opts
        ]
    ]
    def default_hictk_convert_hic_to_mcool_opts = [
        "--tmpdir":"."
    ]
    def hictk_convert_hic_to_mcool_opts = ingest_matrix_opts?.hictk_convert_hic_to_mcool_opts ?: [:]
    def remap = [
        "--resolutions": "-r",
        "--genome": "-g",
        "--verbosity": "-v",
        "--threads": "-t",
        "--compression-lvl": "-l",
        "--force": "-f"
    ]
    def final_hictk_convert_mcool_to_hic_opts = buildCLIOpts(default_hictk_convert_hic_to_mcool_opts, hictk_convert_hic_to_mcool_opts, remap, null)
    logMap += [default_hictk_convert_hic_to_mcool_opts: default_hictk_convert_hic_to_mcool_opts, hictk_convert_hic_to_mcool_opts: hictk_convert_hic_to_mcool_opts]
    def args = [hicFile, output].collect{"'${it}'"}.join(" ")
    def cmd = "hictk convert ${final_hictk_convert_mcool_to_hic_opts} ${args}"
    logMap += [cmd: cmd]
    return [cmd, logMap, output]
}

