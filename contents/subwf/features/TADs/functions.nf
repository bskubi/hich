include {buildCLIOpts} from '../../util/cli.nf'

def buildCmdError(id, mcool, analysisPlan, message) {
    message = "In TADS with id '${id}', mcool '${mcool}', analysisPlan ${analysisPlan}, ${message}"
    error(message)
}

def buildCmd(id, mcool, analysisPlan, cpus) {
    def resolution = analysisPlan?.resolution ?: null
    assert resolution, buildCmdError("analysisPlan.resolution required but was '${resolution}'.")

    def default_hicFindTADs_opts = [
        "--outPrefix": id,
        "--matrix": "${mcool}::/resolutions/${resolution}",
        "-p": cpus,
        "--correctForMultipleTesting": "bonferroni"
    ]
    hicFindTADs_opts = analysisPlan?.hicFindTADs ?: [:]
    def remap = ["--matrix": "-m"]
    hicFindTADs_opts = buildCLIOpts(default_hicFindTADs_opts, hicFindTADs_opts, remap, null)
    def cmd = "hicFindTADs ${hicFindTADs_opts}"

    def output = [
        boundariesBed: "${id}_boundaries.bed",
        boundariesGff: "${id}_boundaries.gff",
        domainsBed: "${id}_domains.bed",
        scoreBedgraph: "${id}_score.bedgraph",
        tadScoreBm: "${id}_tad_score.bm"
    ]

    def logMap = [
        task: "TADS",
        input: [id: id, mcool: mcool, analysisPlan: analysisPlan], 
        output: output
    ]

    return [cmd, logMap, output.boundariesBed, output.boundariesGff, output.domainsBed, output.scoreBedgraph, output.tadScoreBm]
}