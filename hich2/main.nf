// include {ALIGN} from './modules/local/align.nf'
// include {runProcess} from './modules/subworkflows/runProcess.nf'
import hich.Engine

workflow {
    // def sample = [
    //     "condition": "c", 
    //     "biorep": "b", 
    //     "techrep": "t", 
    //     "aligner": "bwa mem", 
    //     "fastq1": "assets/fastq/1k/1k_ERR1413593_1.fq.gz", 
    //     "fastq2": "assets/fastq/1k/1k_ERR1413593_2.fq.gz", 
    //     "aligner_index_dir": "assets/index/bwa/", 
    //     "aligner_index_prefix": "M129"
    // ]

    // channel.of(sample)
    //     | map{Engine.prepSample(SpecKey.ALIGN, it)}
    //     | set{samples}

    // align = runProcess(SpecKey.ALIGN, ALIGN, samples)
    // align.execution_context | map{it.command} | view
}