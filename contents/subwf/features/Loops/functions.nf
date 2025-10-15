include {buildCLIOpts} from '../../util/cli.nf'

def buildCmd(id, matrix, loops_opts) {
    def output = "${id}.loop"
    def logMap = [
        task: "LOOPS", 
        input: [id: id, matrix: matrix, loops_opts: loops_opts],
        output: [loops: output]
    ]

    def default_mustache_opts = [
        "-f": matrix,
        "-o": output,
        "-r": 10000
    ]
    def mustache_opts = loops_opts?.mustache_opts ?: [:]
    logMap += [default_mustache_opts:default_mustache_opts, mustache_opts:mustache_opts]
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
    def final_mustache_opts = buildCLIOpts(default_mustache_opts, mustache_opts, remap, null)
    def cmd = "python -m mustache ${final_mustache_opts}"
    logMap.cmd = cmd

    return [cmd, logMap, output]
}