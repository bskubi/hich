def buildCmd(id, matrix, mustacheParams) {
    def output = "${id}.loop"
    def cmd = ["mustache -f '${matrix}' -o '${output}'"] + mustacheParams
    cmd = cmd.join(" ")
    def logMap = [
        task: "LOOPS", 
        input: [id: id, matrix: matrix, mustacheParams: mustacheParams],
        output: [loops: output]
    ]
    return [cmd, logMap, output]
}