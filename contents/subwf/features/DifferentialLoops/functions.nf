include {buildCLIOpts} from '../../util/cli.nf'

def buildCmdError(params, message) {
    "Error during buildCmd with parameters: ${params}. ${message}"
}

def buildCmd (id1, id2, prefix, matrix1, matrix2, analysisPlan) {
    def output = [
        loop1: "${prefix}.loop1", 
        loop2: "${prefix}.loop2", 
        diffloop1: "${prefix}.diffloop1", 
        diffloop2: "${prefix}.diffloop2"
    ]
    def logMap = [
        task: "DIFFERENTIAL_LOOPS", 
        input: [id1: id1, id2: id2, matrix1: matrix1, matrix2: matrix2, analysisPlan: analysisPlan], 
        output: output
    ]
    def default_diff_mustache_opts = [
        "-f1": matrix1,
        "-f2": matrix2,
        "-o": prefix
    ]
    def diff_mustache_opts = analysisPlan?.diff_mustache ?: [:]
    def remap = [
        "--file1": "-f1",
        "--file2": "-f2",
        "--distance": "-d",
        "--outfile": "-o",
        "--resolution": "-r",
        "--bed1": "-bed1",
        "--matrix1": "-m1",
        "--biases1": "-b1",
        "--bed2": "-bed2",
        "--matrix2": "-m2",
        "--biases2": "-b2",
        "--chromosomeSize": "-cz",
        "--normalization": "-norm",
        "--sparsityThreshold": "-st",
        "--pThreshold": "-pt",
        "--pThreshold2": "-pt2",
        "--sigmaZero": "-sz",
        "--octaves": "-oc",
        "--iterations": "-i",
        "--processes": "-p",
        "--chromosome": "-ch",
        "--chromosome2": "-ch2",
        "--verbose": "-v"
    ]
    diff_mustache_opts = buildCLIOpts(default_diff_mustache_opts, diff_mustache_opts, remap, null)
    def cmd = "diff_mustache ${diff_mustache_opts}"


    return [cmd, logMap, output.loop1, output.loop2, output.diffloop1, output.diffloop2]
}