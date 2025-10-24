import hich.Manifest

process FASTQ_ALIGN {
    container 'bskubi/hich_alignment:latest'
    debug true

    input:
    tuple val(id), 
          path(fastq, arity: '1..2'), 
          val(aligner), 
          path(aligner_index_dir), 
          val(aligner_index_prefix),
          val(config)

    output:
    tuple val(id), path(bam), emit: bam

    script:
    task_plan = Manifest.FastqAlign(id, fastq, aligner, aligner_index_dir, aligner_index_prefix, config, task.cpus)
    bam = task_plan.bam
    task_plan.getScript()

    stub:
    task_plan = Manifest.FastqAlign(id, fastq, aligner, aligner_index_dir, aligner_index_prefix, config, task.cpus)
    bam = task_plan.bam
    task_plan.getStub()
}