package hich
import java.nio.file.Path
import groovy.text.SimpleTemplateEngine
import java.util.Collection

class FastqAlign implements TaskPlan {
    String bam

    FastqAlign ( 
        id, 
        fastq, 
        aligner, 
        aligner_index_dir, 
        aligner_index_prefix, 
        config,
        cpus
    ) {
        this.script_template = '${cmd_align} | ${cmd_samtools_view}'
        this.stub_template = 'touch \'${bam}\''
        this.bam = "${id}.bam"

        fastq = (Collection.isInstance(fastq) ? fastq : [fastq])
                .withIndex()
                .collectEntries{f, i -> ["fastq${i+1}".toString(), f]}

        def bind = [
            id: id, 
            *: fastq, 
            aligner: aligner,
            aligner_index_dir: aligner_index_dir, 
            aligner_index_prefix: aligner_index_prefix,
            cpus: cpus,
            bam: this.bam
        ]

        this.bind_script.cmd_align = new Command(config.cmd_align, bind)
        this.bind_script.cmd_samtools_view = new Command(config.cmd_samtools_view, bind)
        this.bind_stub = [bam: this.bam]
    }

    
}