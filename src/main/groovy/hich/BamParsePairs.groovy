package hich

class BamParsePairs extends TaskPlan {
    String pairs

    BamParsePairs ( 
        id, 
        bam,
        chromsizes,
        config,
        cpus,
        memory
    ) {
        this.script_template = '${command_samtools_view} | ${command_pairtools_parse2} | ${command_pairtools_sort}'
        this.stub_template = 'touch \'${pairs}\''
        this.pairs = "${id}.bam_parse_pairs.pairs.gz"

        def bind = [
            id: id, 
            bam: bam.toString(),
            chromsizes: chromsizes.toString(),
            cpus: cpus,
            memory: memory ? (Math.max(memory.toGiga() - 1, 2)).toString() : "2",
            pairs: this.pairs,
        ]

        this.bind_script.command_samtools_view = new Command(config.command_samtools_view, bind)
        this.bind_script.command_pairtools_parse2 = new Command(config.command_pairtools_parse2, bind)
        this.bind_script.command_pairtools_sort = new Command(config.command_pairtools_sort, bind)
        this.bind_stub = [pairs: this.pairs]
    }
}