include {buildCLIOpts} from '../../util/cli.nf'

def buildCmd(id, pairs, fragmentIndex, tag_restriction_fragments_opts) {
    def logMap = [
        task: "TAG_RESTRICTION_FRAGMENTS",
        input: [
            id: id,
            pairs: pairs,
            fragmentIndex: fragmentIndex,
            tag_restriction_fragments_opts: tag_restriction_fragments_opts
        ]
    ]
    def output = "${id}_fragtag.pairs.gz"
    logMap += [output: [fragPairs: output]]

    def default_hich_pairs_map_ends_opts = [
        "--idx1": "rfrag1",
        "--start1": "rfrag_start1",
        "--end1": "rfrag_end1",
        "--idx2": "rfrag2",
        "--start2": "rfrag_start2",
        "--end2": "rfrag_end2"
    ]
    def hich_pairs_map_ends_opts = tag_restriction_fragments_opts?.hich_pairs_map_ends_opts ?: [:]
    def final_hich_pairs_map_ends_opts = buildCLIOpts(default_hich_pairs_map_ends_opts, hich_pairs_map_ends_opts, [:], null)
    logMap += [default_hich_pairs_map_ends_opts: default_hich_pairs_map_ends_opts, hich_pairs_map_ends_opts: hich_pairs_map_ends_opts]
    def args = [fragmentIndex, pairs, output].collect{"'${it}'"}.join(" ")
    def cmd = "hich pairs map-ends ${final_hich_pairs_map_ends_opts} ${args}"
    logMap.cmd = cmd

    return [cmd, logMap, output]
}