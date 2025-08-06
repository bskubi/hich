include {buildCLIOpts} from '../../util/cli.nf'

def buildCmd(aligner, id, indexDir, indexPrefix, fastq, align_opts, minMapq, cpus) {
    def output = "${id}.bam"
    def inputMap = [
        id: id, 
        fastq: fastq, 
        aligner: aligner, 
        indexDir: indexDir,
        indexPrefix: indexPrefix, 
        align_opts: align_opts,
        minMapq: minMapq
    ]
    def logMap = [task: "ALIGN", output: [bam: output], input: inputMap]
    def alignerCmds = [
        "bwa": "bwa mem",
        "bwa-mem2": "bwa-mem2 mem",
        "bwameth": "bwameth.py",
        "bwameth-mem2": "bwameth.py"
    ]
    def alignerCmd = alignerCmds[aligner] ?: aligner
    def args = fastq.findAll{it}
    def index = "${indexDir}/${indexPrefix}"
    def default_bwa_opts = ["-t": cpus, "-p": true]
    def bwa_opts = null
    def remap = null
    if (aligner in ["bwameth", "bwameth-mem2"]) {
        default_bwa_opts += [
            "--do-not-penalize-chimeras": true,
            "--reference": index
        ]
        bwa_opts = align_opts?.bwameth_opts ?: [:]
        remap = [
            "--threads": "-t",
            "--interleaved": "-p"
        ]
    } else {
        default_bwa_opts += ["-S": true, "-P": true, "-5": true, "-M": true]
        bwa_opts = align_opts?.bwa_mem_opts ?: [:]
        remap = [:]
        args = [index] + args
    }
    args = args.collect{"'${it}'"}
    args = args.join(" ")
    logMap += [default_bwa_opts: default_bwa_opts, bwa_opts: bwa_opts]
    def final_bwa_opts = buildCLIOpts(default_bwa_opts, bwa_opts, remap, null)
    def cmd = "${alignerCmd} ${final_bwa_opts} ${args} | samtools view -b -o '${output}'"

    return [cmd, logMap, output]
}

def getFastq(fastq) {
    return fastq.subMap("fastq", "fastq1", "fastq2").values()
}