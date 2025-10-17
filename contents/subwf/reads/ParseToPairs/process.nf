include {withLog; stubLog} from '../../util/logs.nf'
include {buildCmd} from './functions.nf'
include {validateMemory} from '../../util/memory.nf'

process PARSE_TO_PAIRS {
    publishDir params.general.publish.parse ? params.general.publish.parse : "results",
               saveAs: {params.general.publish.parse ? it : null},
               mode: params.general.publish.mode

    label 'pairs'
    tag "$id"
    conda "$projectDir/env/dev_env.yml"
    container params.general.hichContainer
    debug true

    input:
    tuple val(id), path(sambam), path(chromsizes), val(assembly), val(parse_to_pairs_opts), val(minMapq)

    output:
    tuple val(id), path(output)

    shell:
    memoryV = validateMemory(memory, 2, 2)
    (cmd, logMap, output) = buildCmd(id, sambam, chromsizes, assembly, parse_to_pairs_opts, minMapq, memoryV, task.cpus)
    withLog(cmd, logMap)

    stub:
    memoryV = validateMemory(memory, 2, 2)
    (cmd, logMap, output) = buildCmd(id, sambam, chromsizes, assembly, parse_to_pairs_opts, minMapq, memoryV, task.cpus)
    stub = "touch '${output}'"
    stubLog(stub, cmd, logMap)
}