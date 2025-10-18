include {buildCmd} from './functions.nf'
include {withLog; stubLog} from '../../util/logs.nf'

process ALIGN {
    publishDir params.general.publish.align ? params.general.publish.align : "results",
               saveAs: {params.general.publish.align ? it : null},
               mode: params.general.publish.mode
    
    label 'fullNodeLongTime'
    label 'align'
    tag "$id"
    conda "$projectDir/env/dev_env.yml"
    container params.general.alignmentContainer
    debug true
    
    input:
    tuple val(id), path(fastq), path(fastq1), path(fastq2), val(aligner), path(aligner_index_dir), val(aligner_index_prefix), val(align_opts), val(min_mapq)

    output:
    tuple val(id), path(output)

    shell:
    (cmd, logMap, output) = buildCmd(aligner, id, aligner_index_dir, aligner_index_prefix, fastq, fastq1, fastq2, align_opts, min_mapq, task.cpus)
    withLog(cmd, logMap)

    stub:
    
    (cmd, logMap, output) = buildCmd(aligner, id, aligner_index_dir, aligner_index_prefix, fastq, fastq1, fastq2, align_opts, min_mapq, task.cpus)
    stub = "touch '${output}'"
    stubLog(stub, cmd, logMap)
}