import java.nio.file.Paths
import hich.util.HichUtil
import hich.Engine
import hich.specs.HichSpec.SpecKey
import hich.plans.Plan

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
    val(plan), emit: plan

    script:
    plan = get_plan(
        id, 
        fastq, fastq1, fastq2, 
        aligner, aligner_index_dir, aligner_index_prefix, aligner_opts, 
        task.cpus
    )
    bam = plan.output.bam
    plan.command

    stub:
    plan = get_plan(
        id, 
        fastq, fastq1, fastq2, 
        aligner, aligner_index_dir, aligner_index_prefix, aligner_opts, 
        task.cpus
    )
    bam = context.output.bam
    context.stub
}

Plan get_plan(id, fastq, fastq1, fastq2, aligner, aligner_index_dir, aligner_index_prefix, aligner_opts, cpus) {
    return Engine.getPlan(
        SpecKey.ALIGN,
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