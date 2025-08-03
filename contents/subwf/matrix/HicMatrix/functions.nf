include {validateMemory} from '../../util/memory.nf'

def buildCmd(id, matrixPlanName, pairs, chromsizes, pairsFormat, matrix, juicerToolsPreParams, flags, memory, cpus) {
    juicerToolsPreParams = juicerToolsPreParams ?: []
    matrix = matrix ?: [:]

    // The user can manually set -q in juicerToolsPreParams.
    // Otherwise, use minMapq if specified.
    // Otherwise, no mapq filter is applied.
    if (flags.minMapq instanceof Integer && !juicerToolsPreParams.any{it.contains("-q")}) {
        juicerToolsPreParams += ["-q ${flags.minMapq}"]
    }

    // The user can manually set -r in juicerToolsPreParams.
    // Otherwise, use matrix.resolutions if specified.
    // If neither resolutions nor -r is set, Juicer Tools Pre defaults to producing
    // 2.5M, 1M, 500K, 250K, 100K, 50K, 25K, 10K, and 5K
    if (matrix.resolutions instanceof List && matrix.resolutions && !juicerToolsPreParams.any{it.contains('-r')} ) {
        juicerToolsPreParams += ["-r ${matrix.resolutions.join(',')}"]
    }

    def output = "${id}.hic"
    memory = validateMemory(memory, 2, 8)

    cmd = ["juicer_tools pre -Xms${memory - 2}g -Xmx${memory}g --threads ${cpus}" ] + juicerToolsPreParams + ["'${pairs}' '${output}' '${chromsizes}'"]
    cmd = cmd.findAll{it}
    cmd = cmd.join(" ")
    logMap = [task: "JUICER_TOOLS_PRE", output: output, input: [id: id, pairs: pairs, chromsizes: chromsizes, pairsFormat: pairsFormat, matrix: matrix, juicerToolsPreParams: juicerToolsPreParams, flags: flags]]
    return [cmd, logMap, output]
}