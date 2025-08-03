include {validateMemory} from '../../util/memory.nf'

def buildCmd(id, pairs, chromsizes, sql, cpus, memory) {
    def output = "${id}.pairs.gz"
    def memoryV = validateMemory(memory, 2, 2)
    def reshapeCmd = sql ? ["hich pairs sql --memory-limit '${memoryV}' --threads '${cpus}' '${sql}' '${pairs}'"] : []
    def flipCmd = ["pairtools flip --chroms-path '${chromsizes}'  --nproc-in ${cpus} --nproc-out ${cpus}"]
    
    
    def sortCmd = ["pairtools sort --output '${id}.pairs.gz'  --memory ${memoryV}G --nproc-in ${cpus} --nproc-out ${cpus}"]

    if (!reshapeCmd) {
        flipCmd = flipCmd + ["'${pairs}'"]
        flipCmd = [flipCmd.join(" ")]
    }

    def cmdParts = reshapeCmd + flipCmd + sortCmd
    def cmd = cmdParts.join(" | ")

    def logMap = [
        task: "INGEST_PAIRS", 
        input: [id: id, pairs: pairs, chromsizes: chromsizes, sql: sql], 
        output: [pairs: output]
    ]

    return [cmd, logMap, output]
}