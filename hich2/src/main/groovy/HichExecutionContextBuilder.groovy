import java.nio.file.Path
import groovy.text.SimpleTemplateEngine
import groovy.util.Eval

class HichExecutionContextBuilder {
    static HichExecutionContext alignContext(
        Map commands,
        String id,
        Path fastq,
        Path fastq1,
        Path fastq2,
        String aligner,
        Path aligner_index_dir,
        String aligner_index_prefix,
        Map aligner_opts,
        Integer cpus
    ) {
        HichExecutionContext context = new HichExecutionContext(
            input: [
                id: id,
                fastq: fastq,
                fastq1: fastq1,
                fastq2: fastq2,
                aligner: aligner,
                aligner_index_dir: aligner_index_dir,
                aligner_index_prefix: aligner_index_prefix,
                aligner_opts: aligner_opts,
                cpus: cpus
            ],
            output: [bam:null],
            builder: [
                align: new HichCommandBuilder(),
                samtools_view: new HichCommandBuilder()
            ] 
        )
        def all_errors = []
        if (cpus == null || cpus <= 0 ) {
            all_errors.add("Invalid value of cpus: ${cpus}. Must be positive integer.")
        }

        Map align_CB
        if (["bwa mem", "bwa-mem2"].contains(aligner)) {
            align_CB = commands.bwa_mem_command
        }
        else if (["bwameth", "bwameth-mem2"].contains(aligner)) {
            align_CB = commands.bwameth_command
        }
        else {
            all_errors.add("Unsupported aligner string '${aligner}'")
        }
        
        context.builder.align = (
            context.builder.align
            .setBaseCommand(align_CB.base_command)
            .formatBaseCommand([aligner: aligner])

            .updateOpts(align_CB.default_options)
            .updateArgs(align_CB.arguments)
            .updateOpts(aligner_opts)
            
            
            .evalOpt("-p", [fastq: fastq])
            .formatOpt("-t", [cpus: cpus])
            .formatOpt("--reference", [aligner_index_dir: aligner_index_dir, aligner_index_prefix: aligner_index_prefix])
            .formatArg("aligner_index", [aligner_index_dir: aligner_index_dir, aligner_index_prefix: aligner_index_prefix])
            .evalArg("fastq", [fastq: fastq])
            .evalArg("fastq1", [fastq1: fastq1])
            .evalArg("fastq2", [fastq2: fastq2])
            .buildCommand()
        )

        Map samtools_view_CB = commands.samtools_view_command
        context.builder.samtools_view = (
            context.builder.samtools_view
            .setBaseCommand(samtools_view_CB.base_command)
            .updateOpts(samtools_view_CB.required_options)
            .formatOpt("-o", [id: id])
            .buildCommand()
        )

        String align_command = context.builder.align.prepared_command
        String samtools_view_command = context.builder.samtools_view.prepared_command
        String bam = context.builder.samtools_view.opts["-o"]

        def final_command = "${align_command} | ${samtools_view_command}"
        context.command = final_command
        context.output.bam = bam

        if (all_errors) {
            all_errors.add("FATAL: Errors when building align command.")
            String error_message = all_errors.join("\n")
            throw new Exception(error_message)
        }

        //  Return context to process for execution.
        return context
    }
}