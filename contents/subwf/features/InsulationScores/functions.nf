include {buildCLIOpts} from '../../util/cli.nf'

def buildCmdError(id, mcool, analysisPlan, message) {
    def baseMessage = "In sample with id '${id}', mcool '${mcool}', analysisPlan ${analysisPlan}, ${message}"
    error(baseMessage)
}

def buildCmd(id, mcool, analysisPlan) {
    def cooltools_insulation_opts = analysisPlan?.cooltools_insulation ?: [:]
    def resolution = analysisPlan?.resolution
    def window = analysisPlan?.window
    assert resolution, buildCmdError(id, mcool, analysisPlan, "resolution invalid.")
    assert window, buildCmdError(id, mcool, analysisPlan, "window invalid.")
    
    def tsv = "${id}_insulation.tsv"
    def default_cooltools_insulation_opts = [
        "--bigwig": true,
        "--output": tsv,
    ]
    def remap = [
        "--nproc": "-p",
        "--output": "-o",
        "--regions": "--view"
    ]
    cooltools_insulation_opts = buildCLIOpts(default_cooltools_insulation_opts, cooltools_insulation_opts, remap, null)
    def cool = "${mcool}::/resolutions/${resolution}"
    def cmd = "cooltools insulation ${cooltools_insulation_opts} '${cool}' ${window}"

    def bw = "${tsv}.${window}.bw"
    def output = [
        tsv: tsv, 
        bw: bw
    ]
    def logMap = [
        task: "INSULATION_SCORES", 
        input: [
            id: id, 
            mcool: mcool, 
            analysisPlan: analysisPlan
        ], 
        output: output
    ]
    return [cmd, logMap, tsv, bw]
}