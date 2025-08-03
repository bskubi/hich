include {validateMemory} from '../../util/memory.nf'

def validatePairtoolsParse2Params(parseParams, minMapq) {
    def parseParamsV = parseParams ?: []

    // Use minMapq as default, but override with manually specified --min-mapq
    if (minMapq instanceof Integer && !parseParamsV.any{it.contains("--min-mapq")}) {
        parseParamsV += ["--min-mapq ${minMapq}"]
    } 

    parseParamsV = parseParamsV.join(" ")
}

def buildCmd(id, sambam, chromsizes, assembly, parseParams, sql, minMapq, memory, cpus) {
    def parseParamsV = validatePairtoolsParse2Params(parseParams, minMapq)
    def memoryV = validateMemory(memory, 2, 2)
    def output = "${id}.pairs.gz"
   
    def viewCmd = "samtools view -b '${sambam}'"
    def parse2Cmd = "pairtools parse2 --flip --assembly '${assembly}' --chroms-path '${chromsizes}' ${parseParamsV}"
    def sqlCmd = sql ? "hich pairs sql --memory-limit '${memoryV}' --threads '${cpus}' '${sql}' /dev/stdin" : null
    def pairsSortCmd = "pairtools sort --output '${output}' --memory ${memoryV}G --nproc-in ${cpus} --nproc-out ${cpus}"
    
    def cmd = [viewCmd, parse2Cmd, sqlCmd, pairsSortCmd].findAll{it}.join(" | ")

    def logMap = [
        task: "PARSE",
        input: [
            id: id,
            sambam: sambam,
            chromsizes: chromsizes,
            assembly: assembly,
            parseParams: parseParams,
            sql: sql,
            minMapq: minMapq,
            memory: memory,
            cpus: cpus
        ],
        output: [
            pairs: output
        ]
    ]
    return [cmd, logMap, output]
}
