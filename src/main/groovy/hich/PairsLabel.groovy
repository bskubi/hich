package hich

class PairsLabel extends TaskPlan {
    String pairs_labeled

    PairsLabel ( 
        id, 
        pairs,
        fragment_index,
        config_pairs_label
    ) {
        this.script_template = '${command_hich_pairs_map_ends}'
        this.stub_template = 'touch \'${pairs_labeled}\''
        this.pairs_labeled = "${id}.labeled.pairs.gz"

        def bind = [
            id: id, 
            pairs: pairs.toString(),
            fragment_index: fragment_index.toString(),
            pairs_labeled: this.pairs_labeled
        ]

        this.bind_script.command_hich_pairs_map_ends = new Command(config_pairs_label.command_hich_pairs_map_ends, bind)
        this.bind_stub = [pairs_labeled: this.pairs_labeled]
    }
}