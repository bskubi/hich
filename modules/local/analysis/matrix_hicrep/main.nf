import hich.Hub

process MATRIX_HICREP {
    container 'bskubi/hich_tools:latest'
    
    input:
    tuple val(id), path(matrix), val(command_hich_matrix_hicrep)

    output:
    tuple val(id), path(scc), emit: scc

    script:
    def task_plan = Hub.MatrixHiCRep(id, matrix, command_hich_matrix_hicrep)
    scc = task_plan.scc
    task_plan.getScript()

    stub:
    def task_plan = Hub.MatrixHiCRep(id, matrix, command_hich_matrix_hicrep)
    scc = task_plan.scc
    task_plan.getStub()
}