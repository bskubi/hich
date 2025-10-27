import hich.Hub

process PAIRS_MERGE {
    container 'bskubi/hich_tools:latest'
    
    input:
    tuple val(id), path(pairs), val(config_pairs_merge)

    output:
    tuple val(id), path(pairs), emit: pairs

    script:
    def task_plan = Hub.PairsMerge(id, pairs, config_pairs_merge, task.cpus)
    pairs = task_plan.pairs_merged
    task_plan.getScript()

    stub:
    def task_plan = Hub.PairsMerge(id, pairs, config_pairs_merge, task.cpus)
    pairs = task_plan.pairs_merged
    task_plan.getStub()
}