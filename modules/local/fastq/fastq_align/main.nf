import hich.Hub

process FASTQ_ALIGN {
    container 'bskubi/hich_alignment:latest'
    debug true

    input:
    tuple val(id), 
          path(fastq1),
          path(fastq2),
          path(aligner_index_dir),
          val(config_fastq_align)

    output:
    tuple val(id), path(bam), emit: bam

    script:
    def task_plan = Hub.FastqAlign(id, fastq1, fastq2, aligner_index_dir, config_fastq_align, task.cpus)
    bam = task_plan.bam
    task_plan.getScript()

    stub:
    def task_plan = Hub.FastqAlign(id, fastq1, fastq2, aligner_index_dir, config_fastq_align, task.cpus)
    bam = task_plan.bam
    task_plan.getStub()
}