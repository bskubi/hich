import spock.lang.Specification

// This is the class you are testing.
// Gradle will automatically compile and add your 'src/main/groovy' code to the classpath.
import HichExecutionContextBuilder 
import HichExecutionContext
import static HichContract.PROCESS_INTERFACES
import static HichContract.Schema
import java.nio.file.Paths

class HichExecutionStateBuilderTest extends Specification {

    def "alignContext should return a valid context for BWA alignment of Hi-C data"() {
        def command_blueprint = PROCESS_INTERFACES.ALIGN[Schema.COMMANDS].commands
        def id = "sample1"
        def aligner = "bwa mem"
        def fastq = null
        def fastq1 = Paths.get("R1.fq")
        def fastq2 = Paths.get("R2.fq")
        def aligner_index_dir = Paths.get("index/bwa")
        def aligner_index_prefix = "M129"
        def aligner_opts = ["-T":0]
        def cpus = 1

        when:
        HichExecutionContext context = HichExecutionContextBuilder.alignContext(
            command_blueprint, id, fastq, fastq1, fastq2, aligner, aligner_index_dir, aligner_index_prefix, aligner_opts, cpus
        )

        then:
        context != null
        context.command != null
        context.input != null
        context.output != null
        context.builder != null
        context.builder.align != null
        context.builder.samtools_view != null
        context.input == [
            id: id,
            fastq: fastq,
            fastq1: fastq1,
            fastq2: fastq2,
            aligner: aligner,
            aligner_index_dir: aligner_index_dir,
            aligner_index_prefix: aligner_index_prefix,
            aligner_opts: aligner_opts,
            cpus: cpus
        ]
        context.output == [bam: "${id}.bam"]
        context.builder.align.base_command == "${aligner}"
        context.builder.align.prepared_opts.toSet() == ["-S", "-P", "-5", "-M", "-t '1'", "-T '0'"].toSet()
    }

    def "alignContext should return a valid context for bwameth alignment of interleaved Hi-C data"() {
        def command_blueprint = PROCESS_INTERFACES.ALIGN[Schema.COMMANDS].commands
        def id = "sample1"
        def aligner = "bwameth"
        def fastq = Paths.get("interleaved.fq")
        def fastq1 = null
        def fastq2 = null
        def aligner_index_dir = Paths.get("index/bwameth")
        def aligner_index_prefix = "M129"
        def aligner_opts = ["-T":10]
        def cpus = 5

        when:
        HichExecutionContext context = HichExecutionContextBuilder.alignContext(
            command_blueprint, id, fastq, fastq1, fastq2, aligner, aligner_index_dir, aligner_index_prefix, aligner_opts, cpus
        )

        then:
        context != null
        context.command != null
        context.input != null
        context.output != null
        context.builder != null
        context.builder.align != null
        context.builder.samtools_view != null
        context.input == [
            id: id,
            fastq: fastq,
            fastq1: fastq1,
            fastq2: fastq2,
            aligner: aligner,
            aligner_index_dir: aligner_index_dir,
            aligner_index_prefix: aligner_index_prefix,
            aligner_opts: aligner_opts,
            cpus: cpus
        ]
        context.output == [bam: "${id}.bam"]
        context.builder.align.base_command == "${aligner}"
        def reference = "'${aligner_index_dir}/${aligner_index_prefix}'"
        context.builder.align.prepared_opts.toSet() == ["--do-not-penalize-chimeras", "--reference ${reference}", "-p", "-t '5'", "-T '10'"].toSet()
    }
}