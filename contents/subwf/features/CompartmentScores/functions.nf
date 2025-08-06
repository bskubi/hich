include {buildCLIOpts} from '../../util/cli.nf'

def buildCmd(id, genomeReference, matrix, compartmentScores_opts) {
    def cooltools_eigs_cis_opts = compartmentScores_opts?.cooltools_eigs_cis ?: [:]
    def phasingTrack = "${id}.phase.bed"
    def outPrefix = "${id}_compartments"
    def resolution = compartmentScores_opts?.resolution ?: null
    assert resolution, "Received ${resolution} for required compartmentScores_opts.resolution parameter at which compartment scores are called."

    def default_cooltools_eigs_cis_opts = [
        "--n-eigs": 3,
        "--phasing-track": phasingTrack,
        "--out-prefix": outPrefix,
        "--bigwig": "${matrix}::/resolutions/${resolution}"
    ]
    def remap = [
        "--regions": "--view",
        "--verbose": "-v",
        "--out-prefix": "-o"
    ]
    cooltools_eigs_cis_opts = buildCLIOpts(default_cooltools_eigs_cis_opts, cooltools_eigs_cis_opts, remap, null)
    def cmd = "hich fasta gc ${genomeReference} ${resolution} > ${id}.phase.bed && cooltools eigs-cis ${cooltools_eigs_cis_opts}"
    output = ["${id}.cis.bw", "*.cis.vecs.tsv", "*.cis.lam.txt", "*.phase.bed"]
    
    logMap = [
        task: "COMPARTMENT_SCORES", 
        input: [
            id: id, 
            genomeReference: genomeReference, 
            matrix: matrix, 
            compartmentScores_opts: compartmentScores_opts
        ], 
        output: output
    ]
    return [cmd, logMap, output]
}