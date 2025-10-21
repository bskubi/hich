import java.nio.file.Paths
import hich.util.HichUtil
import hich.Engine

process ALIGN {
    input:
    tuple val (id),
          path(fastq), 
          path(fastq1), 
          path(fastq2),
          val (aligner), 
          path (aligner_index_dir), 
          val (aligner_index_prefix),
          val (aligner_opts)

    output:
    tuple val(id), path(bam), emit: result
    val(execution_context), emit: execution_context

    script:
    context = get_context(
        id, 
        fastq, fastq1, fastq2, 
        aligner, aligner_index_dir, aligner_index_prefix, aligner_opts, 
        task.cpus
    )
    bam = context.output.bam
    execution_context = context
    context.command

    stub:
    context = get_context(
        id, 
        fastq, fastq1, fastq2, 
        aligner, aligner_index_dir, aligner_index_prefix, aligner_opts, 
        task.cpus
    )
    bam = context.output.bam
    execution_context = context
    context.stub
}

def get_context(id, fastq, fastq1, fastq2, aligner, aligner_index_dir, aligner_index_prefix, aligner_opts, cpus) {
    return Engine.getContext(
        hich.specs.HichSpec.ALIGN,
        [
            id, 
            HichUtil.toPath(fastq),
            HichUtil.toPath(fastq1),
            HichUtil.toPath(fastq2),
            aligner,
            HichUtil.toPath(aligner_index_dir),
            aligner_index_prefix,
            aligner_opts,
            cpus
        ]
    )
}