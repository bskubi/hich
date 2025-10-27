package hich

class PairsDedup extends TaskPlan {
    String pairs_deduped

    PairsDedup ( 
        id, 
        pairs,
        config_pairs_dedup,
        cpus
    ) {
        this.script_template = '${command_pairtools_dedup}'
        this.stub_template = 'touch \'${pairs_deduped}\''
        this.pairs_deduped = "${id}.deduped.pairs.gz"

        def bind = [
            id: id, 
            pairs: pairs.toString(),
            pairs_deduped: this.pairs_deduped,
            cpus: cpus
        ]

        this.bind_script.command_pairtools_dedup = new Command(config_pairs_dedup.command_pairtools_dedup, bind)
        this.bind_stub = [pairs_deduped: this.pairs_deduped]
    }
}