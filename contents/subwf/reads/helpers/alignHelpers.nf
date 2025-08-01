include {withLog; stubLog} from '../../util/logs.nf'

def formatFlags(bwaFlags, minMapq) {
    // Convert bwaFlags to list if it's not already (i.e. if user uses "-SP5M" instead of ["-S", "-P", "-5", "-M"])   
    def allFlags = bwaFlags instanceof List ? bwaFlags : [bwaFlags]

    // Add minMapq argument if supplied
    if (minMapq instanceof Integer) {
        allFlags += ["-T ${minMapq}"]
    }

    // Convert non-null arguments to a string and put into single quotes for safety
    def flagsArgs = allFlags.findAll{it}.collect{"'${it}'"}.join(" ")
    return flagsArgs
}

def buildCmdAlignBwa(aligner, id, indexDir, indexPrefix, fastq, bwaFlags, minMapq, cpus) {
    def alignerCmds = [
        "bwa": "bwa mem",
        "bwa-mem2": "bwa-mem2 mem",
        "bwameth": "bwameth.py",
        "bwameth-mem2": "bwameth.py"
    ]
    def flagsArgs = null 
    def alignerCmd = alignerCmds.get(aligner, aligner)
    def fastqArgs = fastq.findAll{it}.collect{"'${it}'"}.join(" ")
    def index = "${indexDir}/${indexPrefix}"
    if (aligner in ["bwameth", "bwameth-mem2"]) {
        flagsArgs = flagsArgs = formatFlags(bwaFlags, null)
        index = "--reference '${index}'"
    } else {
        flagsArgs = formatFlags(bwaFlags, minMapq)
    }
    def cmd = "${alignerCmd} -t ${cpus} ${flagsArgs} ${index} ${fastqArgs} | samtools view -b -o '${id}.bam'"
    def inputMap = [id: id, fastq: fastq, aligner: aligner, index: index, bwaFlags: bwaFlags, minMapq: minMapq]
    def logMap = [task: "AlignBwa", output: "${id}.bam", input: inputMap]
    return [cmd, logMap]
}

def getFastq(fastq) {
    return fastq.subMap("fastq", "fastq1", "fastq2").values()
}