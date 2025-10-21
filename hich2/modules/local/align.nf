import java.nio.file.Paths

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
    HichExecutionContext context = HichEngine.getContext(
        "ALIGN",
        [
            id, 
            HichWorkflowAdapter.toPath(fastq),
            HichWorkflowAdapter.toPath(fastq1),
            HichWorkflowAdapter.toPath(fastq2),
            aligner,
            HichWorkflowAdapter.toPath(aligner_index_dir),
            aligner_index_prefix,
            aligner_opts,
            task.cpus
        ]
    )
    bam = context.output.bam
    execution_context = context
    context.command
}
