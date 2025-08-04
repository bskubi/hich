include {buildCLIOpts} from '../../util/cli.nf'

def buildCmd(aligner, id, indexDir, indexPrefix, fastq, aligner_opts, minMapq, cpus) {
    def alignerCmds = [
        "bwa": "bwa mem",
        "bwa-mem2": "bwa-mem2 mem",
        "bwameth": "bwameth.py",
        "bwameth-mem2": "bwameth.py"
    ]
    def alignerCmd = alignerCmds[aligner] ?: aligner
    def fastqArgs = fastq.findAll{it}.collect{"'${it}'"}.join(" ")
    def index = "${indexDir}/${indexPrefix}"
    def alignerFlags = ["-t": cpus, "-p": true]
    if (aligner in ["bwameth", "bwameth-mem2"]) {
        index = "--reference '${index}'"
        alignerFlags += ["--do-not-penalize-chimeras": true]
        
    } else {
        index = "'${index}'"
        alignerFlags += ["-S": true, "-P": true, "-5": true, "-M": true]
    }
    alignerFlags += (alignerFlags ?: [:])
    def flagsArgs = buildCLIOpts(alignerFlags, null)
    def output = "${id}.bam"
    def cmd = "${alignerCmd} ${flagsArgs} ${index} ${fastqArgs} | samtools view -b -o '${output}'"
    def inputMap = [id: id, fastq: fastq, aligner: aligner, index: index, aligner_opts: aligner_opts, minMapq: minMapq]
    def logMap = [task: "ALIGN", output: "${id}.bam", input: inputMap]
    return [cmd, logMap, output]
}

def getFastq(fastq) {
    return fastq.subMap("fastq", "fastq1", "fastq2").values()
}