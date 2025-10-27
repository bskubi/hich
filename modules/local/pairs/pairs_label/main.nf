import hich.Hub

process PAIRS_LABEL {
    container 'bskubi/hich_tools:latest'
    
    input:
    tuple val(id), path(pairs), path(fragment_index), val(config_pairs_label)

    output:
    tuple val(id), path(pairs_labeled), emit: pairs

    script:
    def task_plan = Hub.PairsLabel(id, pairs, fragment_index, config_pairs_label)
    pairs_labeled = task_plan.pairs_labeled
    task_plan.getScript()

    stub:
    def task_plan = Hub.PairsLabel(id, pairs, fragment_index, config_pairs_label)
    pairs_labeled = task_plan.pairs_labeled
    task_plan.getStub()
}