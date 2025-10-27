package hich

class PairsSelect extends TaskPlan {
    String pairs_selected

    PairsSelect ( 
        id, 
        pairs,
        config_pairs_select,
        cpus
    ) {
        this.script_template = '${command_pairtools_select}'
        this.stub_template = 'touch \'${pairs_selected}\''
        this.pairs_selected = "${id}.selected.pairs.gz"

        def bind = [
            id: id, 
            pairs: pairs.toString(),
            pairs_selected: this.pairs_selected,
            cpus: cpus
        ]

        this.bind_script.command_pairtools_select = new Command(config_pairs_select.command_pairtools_select, bind)
        this.bind_stub = [pairs_selected: this.pairs_selected]
    }
}