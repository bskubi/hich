include {emptyOnLastStep; skip} from '../util/cli.nf'
include {keyJoin} from '../util/keyJoin.nf'
include {withLog; stubLog} from '../util/logs.nf'

process JuicerToolsPre {
    /*
        Juicer Tools Pre documentation: https://github.com/aidenlab/juicer/wiki/Pre

        We use version 1.22.01 as 2.0+ versions are in development and certain
        features available in version 1 are unavailable in 2.
    */
    publishDir params.general.publish.hic ? params.general.publish.hic : "results",
               saveAs: {params.general.publish.hic ? it : null},
               mode: params.general.publish.mode
    container params.general.juicerContainer
    label 'createMatrix'

    input:
    tuple val(id), path(pairs), path(chromsizes), val(pairsFormat), val(matrix), val(juicerToolsPreParams), val(flags)

    output:
    tuple val(id), path(hic)

    shell:
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

    hic = "${id}.hic"
    memory = task.memory ? task.memory.toGiga() : "8"

    cmd = ["java -Xmx${memory}g -jar /app/juicer_tools_1.22.01.jar pre --threads ${task.cpus}" ] + juicerToolsPreParams + ["'${pairs}' '${hic}' '${chromsizes}'"]
    cmd.removeAll([null])
    cmd = cmd.join(" ")
    logMap = [task: "JuicerToolsPre", input: [id: id, pairs: pairs, chromsizes: chromsizes, pairsFormat: pairsFormat, matrix: matrix, juicerToolsPreParams: juicerToolsPreParams, flags: flags], 
    output: [hic: hic]]
    withLog(cmd, logMap)

    stub:
    stub = "touch '${id}.hic'"
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

    hic = "${id}.hic"
    memory = task.memory ? task.memory.toGiga() : "8"

    cmd = ["java -Xmx${memory}g -jar /app/juicer_tools_1.22.01.jar pre --threads ${task.cpus}" ] + juicerToolsPreParams + ["'${pairs}' '${hic}' '${chromsizes}'"]
    cmd.removeAll([null])
    cmd = cmd.join(" ")
    logMap = [task: "JuicerToolsPre", input: [id: id, pairs: pairs, chromsizes: chromsizes, pairsFormat: pairsFormat, matrix: matrix, juicerToolsPreParams: juicerToolsPreParams, flags: flags], 
    output: [hic: hic]]
    stubLog(stub, cmd, logMap)
}

workflow HicMatrix {
    take:
    samples
    
    main:
    samples
        | filter{!skip("hicMatrix") && it.matrix.makeHicFileFormat && (it.pairs || it.latestPairs) && !it.hic}
        | map{tuple(it.id, it.latestPairs, it.chromsizes, it.pairsFormat, it.matrix, it.juicerToolsPreParams, it.subMap("minMapq"))}
        | JuicerToolsPre
        | map{id, hic -> [id: id, hic: hic, latestMatrix: hic]}
        | set{result}
    keyJoin(samples, result, "id") | set{samples}

    samples = emptyOnLastStep("hicMatrix", samples)

    emit:
    samples
}