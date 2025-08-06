include {buildCLIOpts} from '../../util/cli.nf'

def buildCmdError(id, mcool, tads_opts, message) {
    message = "In TADS with id '${id}', mcool '${mcool}', tads_opts ${tads_opts}, ${message}"
    error(message)
}

def buildCmd(id, mcool, tads_opts, cpus) {
    def output = [
        boundariesBed: "${id}_boundaries.bed",
        boundariesGff: "${id}_boundaries.gff",
        domainsBed: "${id}_domains.bed",
        scoreBedgraph: "${id}_score.bedgraph",
        tadScoreBm: "${id}_tad_score.bm"
    ]

    def logMap = [
        task: "TADS",
        input: [id: id, mcool: mcool, tads_opts: tads_opts], 
        output: output
    ]

    def resolution = tads_opts?.resolution ?: null
    assert resolution, buildCmdError(id, mcool, tads_opts, "tads_opts.resolution required but was '${resolution}'.")

    def default_hicfindtads_opts = [
        "--outPrefix": id,
        "--matrix": "${mcool}::/resolutions/${resolution}",
        "-p": cpus,
        "--correctForMultipleTesting": "bonferroni"
    ]
    def hicfindtads_opts = tads_opts?.hicFindTADs ?: [:]
    def remap = ["--matrix": "-m"]
    def final_hicfindtads_opts = buildCLIOpts(default_hicfindtads_opts, hicfindtads_opts, remap, null)
    def cmd = "hicFindTADs ${final_hicfindtads_opts}"



    return [cmd, logMap, output.boundariesBed, output.boundariesGff, output.domainsBed, output.scoreBedgraph, output.tadScoreBm]
}