include {buildCLIOpts} from '../../util/cli.nf'

def buildCmd(id, genomeReference, matrix, compartment_scores_opts) {
    def phasing_track = "${id}.phase.bed"
    output = [bw: "${id}.cis.bw", vecs: "*.cis.vecs.tsv", lam: "*.cis.lam.txt", phasing_track: phasing_track]
    logMap = [
        task: "COMPARTMENT_SCORES", 
        input: [
            id: id, 
            genomeReference: genomeReference, 
            matrix: matrix, 
            compartment_scores_opts: compartment_scores_opts
        ], 
        output: output
    ]

    def out_prefix = "${id}_compartments"
    def resolution = compartment_scores_opts?.resolution ?: null
    assert resolution, "Received '${resolution}' for required compartment_scores_opts.resolution parameter at which compartment scores are called."

    def default_cooltools_eigs_cis_opts = [
        "--n-eigs": 3,
        "--phasing-track": phasing_track,
        "--out-prefix": out_prefix,
        "--bigwig": "${matrix}::/resolutions/${resolution}"
    ]
    def cooltools_eigs_cis_opts = compartment_scores_opts?.cooltools_eigs_cis_opts ?: [:]
    logMap += [default_cooltools_eigs_cis_opts: default_cooltools_eigs_cis_opts, cooltools_eigs_cis_opts: cooltools_eigs_cis_opts]
    def remap = [
        "--regions": "--view",
        "--verbose": "-v",
        "--out-prefix": "-o"
    ]
    def final_cooltools_eigs_cis_opts = buildCLIOpts(default_cooltools_eigs_cis_opts, cooltools_eigs_cis_opts, remap, null)
    def cmd = "hich fasta gc ${genomeReference} ${resolution} > ${phasing_track} && cooltools eigs-cis ${final_cooltools_eigs_cis_opts}"
    logMap.cmd = cmd

    return [cmd, logMap, output]
}