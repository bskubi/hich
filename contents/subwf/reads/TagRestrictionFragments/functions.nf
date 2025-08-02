def buildCmd(id, pairs, fragmentIndex) {
    def output = "${id}_fragtag.pairs.gz"
    def cmd = ["hich pairs map-ends",
           "--idx1 rfrag1 --start1 rfrag_start1 --end1 rfrag_end1",
           "--idx2 rfrag2 --start2 rfrag_start2 --end2 rfrag_end2",
           "'${fragmentIndex}' '${pairs}' '${output}'"
    ]
    cmd = cmd.join(" ")

    def logMap = [
        task: "TAG_RESTRICTION_FRAGMENTS",
        input: [
            id: id,
            pairs: pairs,
            fragmentIndex: fragmentIndex
        ],
        output: [
            pairs: output
        ]
    ]

    return [cmd, logMap, output]
}