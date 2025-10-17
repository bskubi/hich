include {buildCLIOpts} from '../../util/cli.nf'

def _getAlignerConfig(String aligner, String index) {
    def alignerConfigs = [
        "bwa_family": [
            userOptsKey: "bwa_mem_opts",
            defaultOpts: ["-S": true, "-P": true, "-5": true, "-M": true, "-T": 0],
            remap: [:],
            family: "bwa_family"
        ],
        "bwameth_family": [
            userOptsKey: "bwameth_opts",
            defaultOpts: [
                "--do-not-penalize-chimeras": true,
                "--reference": index
            ],
            remap: ["--threads": "-t", "--interleaved": "-p"],
            family: "bwameth_family"
        ]
    ]
    def aligner_families = [
        "bwa": "bwa_family",
        "bwa-mem2": "bwa_family",
        "bwameth": "bwameth_family",
        "bwameth-mem2": "bwameth_family"
    ]
    def family = aligner_families[aligner]
    def config = alignerConfigs[family]
    if (!config) {
        error("Failed to retrieve aligner config for aligner ${aligner} with index ${index}")
    }
    return config
}

def _buildDefaultOpts(Map config, List fastq, minMapq, cpus, Map align_opts) {
    def defaultOpts = config.defaultOpts + ["-t": cpus]

    if (fastq.size() == 1) {
        defaultOpts += ["-p": true]
    }
    if (config.family == "bwa_family" && align_opts?.filterMAPQ && minMapq instanceof Integer) {
        defaultOpts += ["-T": minMapq]
    }
    return defaultOpts
}

def _buildPositionalArgs(Map config, String index, List fastq) {
    if (config.family == "bwa_family") {
        return [index] + fastq
    }
    return fastq // For bwameth, index is an option, not a positional arg
}

def getFastq(fastq) {
    return fastq.subMap("fastq", "fastq1", "fastq2").values()
}

def buildCmd(aligner, id, indexDir, indexPrefix, fastq, align_opts, minMapq, cpus) {
    def output = "${id}.bam"
    def index = "${indexDir}/${indexPrefix}"
    def cleanFastq = fastq.values().findAll { it }

    def alignerCmds = [
        "bwa": "bwa mem",
        "bwa-mem2": "bwa-mem2 mem",
        "bwameth": "python -m bwameth",
        "bwameth-mem2": "python -m bwameth"
    ]
    def alignerCmd = alignerCmds[aligner] ?: aligner

    // 1. Delegate to helper functions
    def config = _getAlignerConfig(aligner, index)
    def defaultOpts = _buildDefaultOpts(config, cleanFastq, minMapq, cpus, align_opts)
    def positionalArgsRaw = _buildPositionalArgs(config, index, cleanFastq)

    // 2. Prepare final command pieces
    def userOpts = align_opts?.get(config.userOptsKey) ?: [:]
    def finalOptsStr = buildCLIOpts(defaultOpts, userOpts, config.remap, null)
    def positionalArgsStr = positionalArgsRaw.collect { "'${it}'" }.join(" ")

    // 3. Assemble the final command and logging map
    def cmd = "${alignerCmd} ${finalOptsStr} ${positionalArgsStr} | samtools view -b -o '${output}'"

    def logMap = [
        task: "ALIGN",
        cmd: cmd,
        output: [bam: output],
        input: [id: id, fastq: fastq, aligner: aligner, index: index],
        bwa_opts_used: [defaults: defaultOpts, user: userOpts, final: finalOptsStr]
    ]

    return [cmd, logMap, output]
}