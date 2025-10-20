import java.nio.file.Path
import groovy.text.SimpleTemplateEngine

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
        Map aligner_opts
    ) {
        
        // User may pass any flags passable to the aligner used.
        aligner_opts = aligner_opts ?: [:]

        def command
        if (["bwa mem", "bwa-mem2"].contains(aligner)) {
            command = commands.bwa_mem_command
        }
        else if (["bwameth", "bwameth-mem2"].containsKey(aligner)) {
            command = commands.bwameth_command
        }
        else {
            throw new Exception("Unsupported aligner string '${aligner}'")
        }
        
        def base_command = command.base_command
        def options = command.default_options + aligner_opts
        def arguments = command.arguments

        def align_base_command = HichUtil.format(base_command, [aligner: aligner])
        
        options = HichUtil.formatIfPresent(options, "-p", [fastq: fastq])
        options["-p"] = options.containsKey("-p") && options["-p"] == "true"
        options = HichUtil.formatIfPresent(options, "-t", [cpus: 5])
        options = HichUtil.formatIfPresent(options, "--reference", [aligner_index_dir: aligner_index_dir, aligner_index_prefix: aligner_index_prefix])

        options = options.findAll{it.value != null && it.value != false}.collect{k, v -> v == true ? k : "${k} '${v}'"}
        
        arguments = HichUtil.formatIfPresent(arguments, "aligner_index", [aligner_index_dir: aligner_index_dir, aligner_index_prefix: aligner_index_prefix])
        arguments = HichUtil.formatIfPresent(arguments, "fastq", [fastq: fastq])
        arguments = HichUtil.formatIfPresent(arguments, "fastq1", [fastq1: fastq1])
        arguments = HichUtil.formatIfPresent(arguments, "fastq2", [fastq2: fastq2])
        arguments = arguments.findAll{it.value != "null"}
        arguments = arguments.findAll{it.value != null}.findAll{it.value != null}.collect{k, v -> "'${v}'"}
        def align_command = ([align_base_command] + options + arguments).join(" ")
        

        command = commands.samtools_view_command
        def samtools_view_base_command = command.base_command
        options = command.required_options
        options = HichUtil.formatIfPresent(options, "-o", [id: id])
        def bam = options["-o"]
        options = options.findAll{it.value != null && it.value != false}.collect{k, v -> v == true ? k : "${k} '${v}'"}
        def samtools_view_command = ([samtools_view_base_command] + options).join(" ")

        def final_command = "${align_command} | ${samtools_view_command}"
        //  Return context to process for execution.
        return new HichExecutionContext(command: final_command, output: [bam: bam])
    }
}