include {validateMemory} from '../../util/memory.nf'
include {buildCLIOpts} from '../../util/cli.nf'

def buildCmd(id, pairs, chromsizes, matrix_opts, minMapq, memory, cpus) {
    def logMap = [
        task: "JUICER_TOOLS_PRE", 
        input: [
            id: id, 
            pairs: pairs, 
            chromsizes: chromsizes, 
            matrix_opts: matrix_opts,
            minMapq: minMapq,
            memory: memory,
            cpus: cpus
        ]
    ]
    def output = "${id}.hic"
    logMap += [output: [hic: output]]
    memory = validateMemory(memory, 2, 8)

    def juicer_tools_pre_opts = matrix_opts?.juicer_tools_pre ?: [:]
    if (minMapq instanceof Integer && !juicer_tools_pre_opts.containsKey("-q")) {
        juicer_tools_pre_opts += ["-q":minMapq]
    }
    def resolutions = matrix_opts?.resolutions ?: [100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000, 500000, 1000000, 2000000, 5000000]
    if (resolutions != null && !juicer_tools_pre_opts["-r"]) {
        resolutions = resolutions instanceof List ? resolutions : [resolutions]
        juicer_tools_pre_opts += ["-r": resolutions.join(',')]
    }

    def default_juicer_tools_pre_opts = [
        "-Xms${memory-2}g": true, 
        "-Xmx${memory}g": true, 
        "--threads": cpus,
        "-t":"."
    ]
    def final_juicer_tools_pre_opts = buildCLIOpts(default_juicer_tools_pre_opts, juicer_tools_pre_opts, [:], null)
    def args = [pairs, output, chromsizes].collect{"'${it}'"}.join(" ")
    def cmd = "juicer_tools pre ${final_juicer_tools_pre_opts} ${args}"

    return [cmd, logMap, output]
}