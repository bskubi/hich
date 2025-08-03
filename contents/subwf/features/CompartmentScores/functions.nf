def buildCmd(id, genomeReference, matrix, resolution, eigsCisParams, nEigs) {
    nEigs = nEigs ?: 3
    eigsCisParams = eigsCisParams ? eigsCisParams.join(" ") : ""
    def cmd = [
        "hich fasta gc ${genomeReference} ${resolution} > ${id}.phase.bed",
        "&& cooltools eigs-cis",
        "--phasing-track ${id}.phase.bed",
        "${eigsCisParams}",
        "--n-eigs ${nEigs}",
        "--out-prefix ${id}_compartments",
        "--bigwig ${matrix}::/resolutions/${resolution}"
    ]
    cmd = cmd.join(" ")

    output = ["${id}.cis.bw", "*.cis.vecs.tsv", "*.cis.lam.txt", "*.phase.bed"]
    
    logMap = [
        task: "COMPARTMENT_SCORES", 
        input: [
            id: id, 
            genomeReference: genomeReference, 
            matrix: matrix, 
            resolution: resolution, 
            eigsCisParams: eigsCisParams
        ], 
        output: output
    ]
    return [cmd, logMap, output]
}