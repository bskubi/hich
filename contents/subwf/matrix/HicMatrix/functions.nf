include {validateMemory} from '../../util/memory.nf'
include {buildCLIOpts} from '../../util/cli.nf'

def buildCmd(id, matrixPlanName, pairs, chromsizes, matrix_opts, minMapq, memory, cpus) {
    def juicer_tools_pre_opts = matrix_opts?.juicer_tools_pre ?: [:]

    if (minMapq instanceof Integer && !juicer_tools_pre_opts.contains("-q")) {
        juicer_tools_pre_opts += ["-q":minMapq]
    }

    def resolutions = matrix_opts?.resolutions
    if (resolutions != null && !juicer_tools_pre_opts["-r"]) {
        resolutions = resolutions instanceof List ? resolutions : [resolutions]
        juicer_tools_pre_opts += ["-r": resolutions.join(',')]
    }

    def output = "${id}.hic"
    memory = validateMemory(memory, 2, 8)

    def juicer_tools_pre_defaults = [
        "-Xms${memory-2}g": true, 
        "-Xmx${memory}g": true, 
        "--threads": cpus
    ]
    juicer_tools_pre_opts = buildCLIOpts(juicer_tools_pre_defaults, juicer_tools_pre_opts)
    def cmd = "juicer_tools pre ${juicer_tools_pre_opts} '${pairs}' '${output}' '${chromsizes}'"
    def logMap = [
        task: "JUICER_TOOLS_PRE", 
        output: [hic: output], 
        input: [
            id: id, 
            pairs: pairs, 
            chromsizes: chromsizes, 
            matrix_opts: matrix_opts,
            minMapq: minMapq
        ]
    ]
    return [cmd, logMap, output]
}