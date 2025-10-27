package hich

class PairsMerge extends TaskPlan {
    String pairs_merged

    PairsMerge ( 
        id, 
        pairs,
        config_pairs_merge,
        cpus
    ) {
        this.script_template = '${command_pairtools_merge}'
        this.stub_template = 'touch \'${pairs_merged}\''
        this.pairs_merged = "${id}.merged.pairs.gz"

        def bind = [
            id: id, 
            pairs: pairs.collect{it.toString()}.join(" "),
            pairs_merged: this.pairs_merged,
            cpus: cpus
        ]

        this.bind_script.command_pairtools_merge = new Command(config_pairs_merge.command_pairtools_merge, bind)
        this.bind_stub = [pairs_merged: this.pairs_merged]
    }
}