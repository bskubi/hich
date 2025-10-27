import hich.Hub

process PAIRS_DEDUP {
    container 'bskubi/hich_tools:latest'
    
    input:
    tuple val(id), path(pairs), val(config_pairs_dedup)

    output:
    tuple val(id), path(pairs), emit: pairs

    script:
    def task_plan = Hub.PairsDedup(id, pairs, config_pairs_dedup, task.cpus)
    pairs = task_plan.pairs_deduped
    task_plan.getScript()

    stub:
    def task_plan = Hub.PairsDedup(id, pairs, config_pairs_dedup, task.cpus)
    pairs = task_plan.pairs_deduped
    task_plan.getStub()
}