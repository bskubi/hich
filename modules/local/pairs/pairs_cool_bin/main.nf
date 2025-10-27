import hich.Hub
process PAIRS_COOL_BIN {
    container 'bskubi/hich_cooler:latest'

    input:
    tuple val(id), path(pairs), path(chromsizes), val(config_pairs_cool_bin)

    output:
    tuple val(id), path(cool), emit: cool

    script:
    def task_plan = Hub.PairsCoolBin(id, pairs, chromsizes, config_pairs_cool_bin)
    cool = task_plan.cool
    task_plan.getScript()
    
    stub:
    def task_plan = Hub.PairsCoolBin(id, pairs, chromsizes, config_pairs_cool_bin)
    cool = task_plan.cool
    task_plan.getStub()
}