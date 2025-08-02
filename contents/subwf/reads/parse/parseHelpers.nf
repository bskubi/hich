def validateParseParams(parseParams, flags) {
    def parseParamsV = parseParams ?: []

    // Use minMapq as default, but override with manually specified --min-mapq
    if (flags.minMapq instanceof Integer && !parseParamsV.any{it.contains("--min-mapq")}) {
        parseParamsV += ["--min-mapq ${flags.minMapq}"]
    } 

    parseParamsV = parseParamsV.join(" ")
}

def validateMemory(memory) {
    // Use 2G less than the total memory allocated for the job
    // to a minimum of 2G
    memory ? Math.max(memory.toGiga() - 2, 2) : 2
}

def buildCmdPairtoolsParse2(id, sambam, chromsizes, assembly, parseParams, sql, flags, memory, cpus) {
    def parseParamsV = validateParseParams(parseParams, flags)
    def memoryV = validateMemory(memory)
   
    def viewCmd = "samtools view -b '${sambam}'"
    def parse2Cmd = "pairtools parse2 --flip --assembly '${assembly}' --chroms-path '${chromsizes}' ${parseParamsV}"
    def sqlCmd = sql ? "hich pairs sql --memory-limit '${memoryV}' --threads '${cpus}' '${sql}' /dev/stdin" : null
    def pairsSortCmd = "pairtools sort --output '${id}.pairs.gz' --memory ${memoryV}G --nproc-in ${cpus} --nproc-out ${cpus}"
    
    def cmd = [viewCmd, parse2Cmd, sqlCmd, pairsSortCmd].findAll{it}.join(" | ")

    def logMap = [
        task: "PairtoolsParse2",
        input: [
            id: id,
            sambam: sambam,
            chromsizes: chromsizes,
            assembly: assembly,
            parseParams: parseParams,
            sql: sql,
            flags: flags,
            memory: memory,
            cpus: cpus
        ],
        output: [
            pairs: "${id}.pairs.gz"
        ]
    ]
    return [cmd, logMap]
}
