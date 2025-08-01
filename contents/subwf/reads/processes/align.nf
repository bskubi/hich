include {buildCmdAlignBwa} from '../helpers/alignHelpers.nf'
include {withLog; stubLog} from '../../util/logs.nf'

process Align {
    publishDir params.general.publish.align ? params.general.publish.align : "results",
               saveAs: {params.general.publish.align ? it : null},
               mode: params.general.publish.mode
    
    label 'whenLocal_allConsuming'
    label 'align'
    tag "$id"
    conda "$projectDir/env/dev_env.yml"
    container params.general.alignmentContainer
    
    input:
    tuple val(id), val(aligner), path(indexDir), val(indexPrefix), path(fastq, arity: 1..2), val(bwaFlags), val(minMapq)

    output:
    tuple val(id), path("${id}.bam")

    shell:
    (cmd, logMap) = buildCmdAlignBwa(aligner, id, indexDir, indexPrefix, fastq, bwaFlags, minMapq, task.cpus)
    withLog(cmd, logMap)

    stub:
    stub = "touch '${id}.bam'"
    (cmd, logMap) = buildCmdAlignBwa(aligner, id, indexDir, indexPrefix, fastq, bwaFlags, minMapq, task.cpus)
    stubLog(stub, cmd, logMap)
}