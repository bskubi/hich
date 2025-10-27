package hich

class PairsCoolBin extends TaskPlan {
    String cool

    PairsCoolBin ( 
        id, 
        pairs,
        chromsizes, 
        config_pairs_cool_bin
    ) {
        this.script_template = '${command_cooler_cload_pairs}'
        this.stub_template = 'touch \'${cool}\''
        this.cool = "${id}.cool"

        def bind = [
            id: id, 
            pairs: pairs.toString(),
            chromsizes: chromsizes.toString(),
            cool: this.cool
        ]

        this.bind_script.command_cooler_cload_pairs = new Command(config_pairs_cool_bin.command_cooler_cload_pairs, bind)
        this.bind_stub = [cool: this.cool]
    }
}