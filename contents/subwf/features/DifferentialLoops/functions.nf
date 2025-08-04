def buildCmd (id1, id2, prefix, matrix1, matrix2, mustacheParams) {
    def cmd = [
        "diff_mustache",
        "-f1 '${matrix1}'",
        "-f2 '${matrix2}'",
        "-o '${prefix}'"
    ] + mustacheParams
    cmd = cmd.findAll{it}
    cmd = cmd.join(" ")
    def output = [
        loop1: "${prefix}.loop1", 
        loop2: "${prefix}.loop2", 
        diffloop1: "${prefix}.diffloop1", 
        diffloop2: "${prefix}.diffloop2"
    ]
    def logMap = [
        task: "DIFFERENTIAL_LOOPS", 
        input: [id1: id1, id2: id2, matrix1: matrix1, matrix2: matrix2, mustacheParams: mustacheParams], 
        output: output
    ]
    return [cmd, logMap, output.loop1, output.loop2, output.diffloop1, output.diffloop2]
}