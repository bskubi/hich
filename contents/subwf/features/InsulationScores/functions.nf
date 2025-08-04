def buildCmd(id, mcool, resolution, cooltoolsInsulationParams, window) {
    def tsv = "${id}_insulation.tsv"
    def bw = "${tsv}.${window}.bw"
    def cmd = [
        "cooltools insulation",
        "--bigwig",
        "--output '${tsv}'"] +
          cooltoolsInsulationParams +
          ["'${mcool}::/resolutions/${resolution}' ${window}"]
    cmd = cmd.findAll{it}
    
    cmd = cmd.join(" ")
    def output = [
        tsv: tsv, 
        bw: bw
    ]
    def logMap = [
        task: "INSULATION_SCORES", 
        input: [
            id: id, 
            mcool: mcool, 
            resolution: resolution, 
            cooltoolsInsulationParams: cooltoolsInsulationParams, 
            window: window
        ], 
        output: output
    ]
    return [cmd, logMap, tsv, bw]
}