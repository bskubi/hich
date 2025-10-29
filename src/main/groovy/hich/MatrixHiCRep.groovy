package hich

class MatrixHiCRep extends TaskPlan {
    String scc

    MatrixHiCRep ( 
        id, 
        matrix,
        command_hich_matrix_hicrep
    ) {
        this.script_template = '${command_hich_matrix_hicrep}'
        this.stub_template = 'touch \'${scc}\''
        this.scc = "${id}.scc.txt"

        def bind = [
            id: id,
            scc: this.scc,
            matrix: matrix.join(" ")
        ]

        this.bind_script.command_hich_matrix_hicrep = new Command(command_hich_matrix_hicrep, bind)
        this.bind_stub = [scc: this.scc]
    }
}