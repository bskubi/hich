package hich

class FastqAlign extends TaskPlan {
    String bam

    FastqAlign ( 
        id, 
        fastq1,
        fastq2,
        aligner_index_dir,
        config,
        cpus
    ) {
        this.script_template = '${command_align} | ${command_samtools_view}'
        this.stub_template = 'touch \'${bam}\''
        this.bam = "${id}.bam"

        def bind = [
            id: id, 
            fastq1: fastq1.toString(),
            fastq2: fastq2.toString(),
            aligner_index_dir: aligner_index_dir.toString(),
            cpus: cpus,
            bam: this.bam
        ]

        this.bind_script.command_align = new Command(config.command_align, bind)
        this.bind_script.command_samtools_view = new Command(config.command_samtools_view, bind)
        this.bind_stub = [bam: this.bam]
    }
}