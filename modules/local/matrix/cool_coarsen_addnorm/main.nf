import hich.Hub

process COOL_COARSEN_ADDNORM {
    container 'bskubi/hich_cooler:latest'

    input:
    tuple val(id), path(cool), val(config_cool_coarsen_addnorm)

    output:
    tuple val(id), path(mcool), emit: mcool

    script:
    def task_plan = Hub.CoolCoarsenAddNorm(id, cool, config_cool_coarsen_addnorm, task.cpus)
    mcool = task_plan.mcool
    task_plan.getScript()

    stub:
    def task_plan = Hub.CoolCoarsenAddNorm(id, cool, config_cool_coarsen_addnorm, task.cpus)
    mcool = task_plan.mcool
    task_plan.getStub()
}