def buildCmdMcoolToHic(id, mcoolFile, cpus) {
    def output = "${id}.hic"
    cpus = cpus >= 2 ? cpus : 2
    threads = "-t ${cpus}"
    cmd = "mkdir ./tmp && hictk convert --tmpdir ./tmp ${threads} '${mcoolFile}' '${output}'"
    def logMap = [task: "MCOOL_TO_HIC", output: [hic: output], input: [id: id, mcoolFile: mcoolFile]]
    return [cmd, logMap, output]
}

def buildCmdHicToMcool(id, hicFile) {
    def output = "${id}.mcool"
    def cmd = "mkdir ./tmp && hictk convert --tmpdir ./tmp '${hicFile}' '${output}'"
    def logMap = [task: "HIC_TO_MCOOL", output: [mcool: output], input: [id: id, hicFile: hicFile]]
    return [cmd, logMap, output]
}

