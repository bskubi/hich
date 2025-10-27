import hich.Hub

process PAIRS_SELECT {
    container 'bskubi/hich_tools:latest'
    
    input:
    tuple val(id), path(pairs), val(config_pairs_select)

    output:
    tuple val(id), path(pairs), emit: pairs

    script:
    def task_plan = Hub.PairsSelect(id, pairs, config_pairs_select, task.cpus)
    pairs = task_plan.pairs_selected
    task_plan.getScript()


    stub:
    def task_plan = Hub.PairsSelect(id, pairs, config_pairs_select, task.cpus)
    pairs = task_plan.pairs_selected
    task_plan.getStub()
}