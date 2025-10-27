import hich.Hub

process BAM_PARSE_PAIRS {
    container 'bskubi/hich_tools:latest'
    memory 4.GB
    
    input:
    tuple val(id), path(bam), path(chromsizes), val(config_bam_parse_pairs)

    output:
    tuple val(id), path(pairs), emit: pairs

    script:
    def task_plan = Hub.BamParsePairs(id, bam, chromsizes, config_bam_parse_pairs, task.cpus, task.memory)
    pairs = task_plan.pairs
    task_plan.getScript()

    stub:
    def task_plan = Hub.BamParsePairs(id, bam, chromsizes, config_bam_parse_pairs, task.cpus, task.memory)
    pairs = task_plan.pairs
    task_plan.getStub()
}