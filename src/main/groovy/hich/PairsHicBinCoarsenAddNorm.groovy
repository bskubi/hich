package hich

class PairsHicBinCoarsenAddNorm extends TaskPlan {
    String hic

    PairsHicBinCoarsenAddNorm ( 
        id, 
        pairs,
        chromsizes, 
        config_pairs_hic_bin_coarsen_addnorm,
        cpus,
        memory
    ) {
        this.script_template = '${command_juicer_tools_pre}'
        this.stub_template = 'touch \'${hic}\''
        this.hic = "${id}.hic"

        def bind = [
            id: id, 
            pairs: pairs.toString(),
            chromsizes: chromsizes.toString(),
            hic: this.hic,
            cpus: cpus,
            xmx: memory ? (Math.max(memory.toGiga() - 1, 2)).toString() : "2",
            xms: memory ? (Math.max(memory.toGiga() - 2, 1)).toString() : "1",
        ]

        this.bind_script.command_juicer_tools_pre = new Command(config_pairs_hic_bin_coarsen_addnorm.command_juicer_tools_pre, bind)
        this.bind_stub = [hic: this.hic]
    }
}