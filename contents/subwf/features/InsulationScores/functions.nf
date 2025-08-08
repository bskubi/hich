include {buildCLIOpts} from '../../util/cli.nf'

def buildCmdError(id, mcool, insulation_scores_opts, message) {
    def baseMessage = "In sample with id '${id}', mcool '${mcool}', insulation_scores_opts ${insulation_scores_opts}: ${message}"
    error(baseMessage)
}

def buildCmd(id, mcool, insulation_scores_opts) {
    def tsv = "${id}_insulation.tsv"
    def resolution = insulation_scores_opts?.resolution ?: 10000
    def window = insulation_scores_opts?.window ?: resolution * 10
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
            insulation_scores_opts: insulation_scores_opts
        ], 
        output: output
    ]
    def cooltools_insulation_opts = insulation_scores_opts?.cooltools_insulation_opts ?: [:]
    
    
    assert resolution, buildCmdError(id, mcool, insulation_scores_opts, "resolution invalid.")
    assert window, buildCmdError(id, mcool, insulation_scores_opts, "window invalid.")
    def default_cooltools_insulation_opts = [
        "--bigwig": true,
        "--output": tsv,
    ]
    def remap = [
        "--nproc": "-p",
        "--output": "-o",
        "--regions": "--view"
    ]
    logMap += [default_cooltools_insulation_opts:default_cooltools_insulation_opts, cooltools_insulation_opts:cooltools_insulation_opts]
    def final_cooltools_insulation_opts = buildCLIOpts(default_cooltools_insulation_opts, cooltools_insulation_opts, remap, null)
    def cool = "${mcool}::/resolutions/${resolution}"
    def args = [cool, window].collect{"'${it}'"}.join(" ")
    def cmd = "cooltools insulation ${final_cooltools_insulation_opts} ${args}"

    return [cmd, logMap, tsv, bw]
}