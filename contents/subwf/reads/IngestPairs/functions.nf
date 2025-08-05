include {validateMemory} from '../../util/memory.nf'
include {buildCLIOpts} from '../../util/cli.nf'

def buildCmd(id, pairs, chromsizes, ingestPairs_opts, cpus, memory) {
    def memoryV = validateMemory(memory, 2, 2)
    
    def hichPairsReshape_opts = ["--memory-limit":memoryV, "--threads":cpus]
    def sql = ingestPairs_opts?.hich_pairs_reshape?.sql
    hichPairsReshape_opts = buildCLIOpts(hichPairsReshape_opts, ingestPairs_opts?.hich_pairs_reshape)
    def reshapeCmd = sql ? "hich pairs sql ${hichPairsReshape_opts} '${sql}' '${pairs}'" : null
    
    def pairtoolsFlip_opts = ["--chroms-path":chromsizes, "--nproc-in":cpus, "--nproc-out":cpus]
    pairtoolsFlip_opts = buildCLIOpts(pairtoolsFlip_opts, ingestPairs_opts?.pairtools_flip)
    def flipCmd = "pairtools flip ${pairtoolsFlip_opts}" + (reshapeCmd ? "" : " '${pairs}'")

    def output = "${id}.pairs.gz"
    def pairtoolsSort_opts = ["--output":output, "--memory":"${memoryV}G", "--nproc-in":cpus, "--nproc-out":cpus]
    pairtoolsSort_opts = buildCLIOpts(pairtoolsSort_opts, ingestPairs_opts?.pairtools_sort)
    def sortCmd = "pairtools sort ${pairtoolsSort_opts}"

    def cmdParts = [reshapeCmd, flipCmd, sortCmd]
    cmdParts = cmdParts.findAll{it}
    def cmd = cmdParts.join(" | ")

    def logMap = [
        task: "INGEST_PAIRS", 
        input: [id: id, pairs: pairs, chromsizes: chromsizes, ingestPairs_opts: ingestPairs_opts], 
        output: [pairs: output]
    ]

    return [cmd, logMap, output]
}