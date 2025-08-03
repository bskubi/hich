def buildCmd(id, mcool, resolution, hicExplorerFindTADsParams, cpus) {
    def cmd = ["hicFindTADs --outPrefix ${id} --matrix ${mcool}::/resolutions/${resolution} -p ${cpus}"] + hicExplorerFindTADsParams
    cmd = cmd.join(" ")

    def output = [
        boundariesBed: "${id}_boundaries.bed",
        boundariesGff: "${id}_boundaries.gff",
        domainsBed: "${id}_domains.bed",
        scoreBedgraph: "${id}_score.bedgraph",
        tadScoreBm: "${id}_tad_score.bm"
    ]

    def logMap = [
        task: "TADS",
        input: [id: id, mcool: mcool, resolution: resolution, hicExplorerFindTADsParams: hicExplorerFindTADsParams], 
        output: output
    ]

    return [cmd, logMap, output.boundariesBed, output.boundariesGff, output.domainsBed, output.scoreBedgraph, output.tadScoreBm]
}