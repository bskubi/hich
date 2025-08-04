include {buildCmd} from './functions.nf'
include {withLog; stubLog} from '../../util/logs.nf'

process ALIGN {
    publishDir params.general.publish.align ? params.general.publish.align : "results",
               saveAs: {params.general.publish.align ? it : null},
               mode: params.general.publish.mode
    
    label 'whenLocal_allConsuming'
    label 'align'
    tag "$id"
    conda "$projectDir/env/dev_env.yml"
    container params.general.alignmentContainer
    debug true
    
    input:
    tuple val(id), val(aligner), path(indexDir), val(indexPrefix), path(fastq, arity: 1..2), val(bwaFlags), val(minMapq)

    output:
    tuple val(id), path(output)

    shell:
    (cmd, logMap, output) = buildCmd(aligner, id, indexDir, indexPrefix, fastq, bwaFlags, minMapq, task.cpus)
    withLog(cmd, logMap)

    stub:
    
    (cmd, logMap, output) = buildCmd(aligner, id, indexDir, indexPrefix, fastq, bwaFlags, minMapq, task.cpus)
    stub = "touch '${output}'"
    stubLog(stub, cmd, logMap)
}