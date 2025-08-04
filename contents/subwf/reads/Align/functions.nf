include {buildFlags} from '../../util/cli.nf'

def buildCmd(aligner, id, indexDir, indexPrefix, fastq, bwaFlags, minMapq, cpus) {
    def alignerCmds = [
        "bwa": "bwa mem",
        "bwa-mem2": "bwa-mem2 mem",
        "bwameth": "bwameth.py",
        "bwameth-mem2": "bwameth.py"
    ]
    def alignerCmd = alignerCmds.get(aligner, aligner)
    def fastqArgs = fastq.findAll{it}.collect{"'${it}'"}.join(" ")
    def index = "${indexDir}/${indexPrefix}"
    bwaFlags += ["-t": cpus]
    if (aligner in ["bwameth", "bwameth-mem2"]) {
        index = "--reference '${index}'"
    } else {
        bwaFlags += ["-t": minMapq]
        index = "'${index}'"
    }
    def flagsArgs = buildFlags(bwaFlags)
    def output = "${id}.bam"
    def cmd = "${alignerCmd} ${flagsArgs} ${index} ${fastqArgs} | samtools view -b -o '${output}'"
    def inputMap = [id: id, fastq: fastq, aligner: aligner, index: index, bwaFlags: bwaFlags, minMapq: minMapq]
    def logMap = [task: "ALIGN", output: "${id}.bam", input: inputMap]
    return [cmd, logMap, output]
}

def getFastq(fastq) {
    return fastq.subMap("fastq", "fastq1", "fastq2").values()
}