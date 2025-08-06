include {validateMemory} from '../../util/memory.nf'
include {buildCLIOpts} from '../../util/cli.nf'

def buildCmd(id, pairs, chromsizes, ingest_pairs_opts, cpus, memory) {
    def logMap = [
        task: "INGEST_PAIRS", 
        input: [
            id: id, 
            pairs: pairs, 
            chromsizes: chromsizes, 
            ingest_pairs_opts: ingest_pairs_opts,
            cpus: cpus,
            memory: memory
        ]
    ]
    def memoryV = validateMemory(memory, 2, 2)
    
    def default_hich_pairs_sql_opts = ["--memory-limit":memoryV, "--threads":cpus]
    def hich_pairs_sql_opts = ingest_pairs_opts?.hich_pairs_reshape ?: [:]
    def sql = hich_pairs_sql_opts?.sql
    def final_hich_pairs_sql_opts = buildCLIOpts(default_hich_pairs_sql_opts, hich_pairs_sql_opts, [:], null)
    logMap += [default_hich_pairs_sql_opts: default_hich_pairs_sql_opts, hich_pairs_sql_opts: hich_pairs_sql_opts]
    def hich_pairs_sql_cmd = sql ? "hich pairs sql ${final_hich_pairs_sql_opts} '${sql}' '${pairs}'" : null
    
    def default_pairtools_flip_opts = ["--chroms-path":chromsizes, "--nproc-in":cpus, "--nproc-out":cpus]
    def pairtools_flip_opts = ingest_pairs_opts?.pairtools_flip_opts ?: [:]
    def remap = [
        "--chroms-path": "-c",
        "--output": "-o"
    ]
    final_pairtools_flip_opts = buildCLIOpts(default_pairtools_flip_opts, pairtools_flip_opts, remap, null)
    logMap += [default_pairtools_flip_opts: default_pairtools_flip_opts, pairtools_flip_opts: pairtools_flip_opts]
    def pairtools_flip_cmd = "pairtools flip ${final_pairtools_flip_opts}" + (hich_pairs_sql_cmd ? "" : " '${pairs}'")

    def output = "${id}.pairs.gz"
    logMap += [output: [pairs: output]]
    def default_pairtools_sort_opts = ["--output":output, "--memory":"${memoryV}G", "--nproc-in":cpus, "--nproc-out":cpus]
    def pairtools_sort_opts = ingest_pairs_opts?.pairtools_sort_opts ?: [:]
    final_pairtools_sort_opts = buildCLIOpts(default_pairtools_sort_opts, pairtools_sort_opts, [:], null)
    def pairtools_sort_cmd = "pairtools sort ${final_pairtools_sort_opts}"

    def cmd = [hich_pairs_sql_cmd, pairtools_flip_cmd, pairtools_sort_cmd]
    cmd = cmd.findAll{it}
    cmd = cmd.join(" | ")
    logMap.cmd = cmd

    return [cmd, logMap, output]
}