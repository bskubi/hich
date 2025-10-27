import hich.Hub

process PAIRS_HIC_BIN_COARSEN_ADDNORM {
    container 'bskubi/hich_juicer:latest'

    input:
    tuple val(id), path(pairs), path(chromsizes), val(config_pairs_hic_bin_coarsen_addnorm)

    output:
    tuple val(id), path(hic), emit: hic

    script:
    def task_plan = Hub.PairsHicBinCoarsenAddNorm(id, pairs, chromsizes, config_pairs_hic_bin_coarsen_addnorm, task.cpus, task.memory)
    hic = task_plan.hic
    task_plan.getScript()

    stub:
    def task_plan = Hub.PairsHicBinCoarsenAddNorm(id, pairs, chromsizes, config_pairs_hic_bin_coarsen_addnorm, task.cpus, task.memory)
    hic = task_plan.hic
    task_plan.getStub()
    
}