include {buildCLIOpts} from '../../util/cli.nf'

def buildCmd(id, matrix, analysisPlan) {
    def output = "${id}.loop"
    def default_mustache_opts = [
        "-f": matrix,
        "-o": output
    ]
    def remap = [
        "--file": "-f",
        "--distance": "-d",
        "--outfile": "-o",
        "--resolution": "-r",
        "--bed": "-bed",
        "--matrix": "-m",
        "--biases": "-b",
        "--chromosomeSize": "-cz",
        "--normalization": "-norm",
        "--sparsityThreshold": "-st",
        "--pThreshold": "-pt",
        "--sigmaZero": "-sz",
        "--octaves": "-oc",
        "--processes": "-p",
        "--chromosome": "-ch",
        "--chromosome2": "-ch2",
        "--verbose": "-v"
    ]
    def mustache_opts = analysisPlan?.mustache ?: [:]
    mustache_opts = buildCLIOpts(default_mustache_opts, mustache_opts, remap, null)
    def cmd = "mustache ${mustache_opts}"
    def logMap = [
        task: "LOOPS", 
        input: [id: id, matrix: matrix, analysisPlan: analysisPlan],
        output: [loops: output]
    ]
    return [cmd, logMap, output]
}