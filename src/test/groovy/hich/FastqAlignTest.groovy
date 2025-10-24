import spock.lang.Specification
import hich.FastqAlign
import groovy.yaml.YamlSlurper
import java.nio.file.Paths

class FastqAlignSpec extends Specification {

    // Define the class under test
    FastqAlign fastq_align 

    // A blank feature method (test)
    def "Test stub command"() {
        setup:
        def hich_config_file = Paths.get("hich_config.yaml")
        
        when:
        def hich_config = new YamlSlurper().parse(hich_config_file)
        def fastq_align = new FastqAlign(1, ["1.fq", "2.fq"], "bwa", "tests/assets/bwa", "M129", hich_config.standard_alignment.config_fastq_align, 5)
        
        then:
        hich_config.standard_alignment.config_fastq_align.cmd_align != null
        hich_config.standard_alignment.config_fastq_align.cmd_samtools_view != null
        fastq_align.getStub() == "touch '1.bam'"
    }
}