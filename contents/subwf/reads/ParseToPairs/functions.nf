include {validateMemory} from '../../util/memory.nf'
include {buildCLIOpts} from '../../util/cli.nf'


def buildCmd(id, sambam, chromsizes, assembly, parse_to_pairs_opts, minMapq, memory, cpus) {
    def logMap = [
        task: "PARSE_TO_PAIRS",
        input: [
            id: id,
            sambam: sambam,
            chromsizes: chromsizes,
            assembly: assembly,
            parse_to_pairs_opts: parse_to_pairs_opts,
            minMapq: minMapq,
            memory: memory,
            cpus: cpus
        ]
    ]

    def memoryV = validateMemory(memory, 2, 2)
    def output = "${id}.pairs.gz" 
    logMap += [output: [pairs: output]]

    def samtoolsViewCmd = "samtools view -b '${sambam}'"

    def default_pairtools_parse2_opts = [
        "--assembly": assembly, 
        "--chroms-path": chromsizes, 
        "--min-mapq": minMapq, 
        "--flip": true
    ]
    def pairtools_parse2_opts = parse_to_pairs_opts?.pairtools_parse2_opts ?: [:]
    logMap += [default_pairtools_parse2_opts: default_pairtools_parse2_opts, pairtools_parse2_opts: pairtools_parse2_opts]
    def remap = ["--chroms-path": "-c"]
    def final_pairtools_parse2_opts = buildCLIOpts(default_pairtools_parse2_opts, pairtools_parse2_opts, remap, null)
    def pairtools_parse2_cmd = "pairtools parse2 ${final_pairtools_parse2_opts}"
    
    def hich_pairs_sql_cmd = null
    if (parse_to_pairs_opts?.hichPairsSql?.sql) {
        def default_hich_pairs_sql_opts = [
            "--memory-limit": memoryV, 
            "--threads": cpus
        ]
        def hich_pairs_sql_opts = parse_to_pairs_opts?.hich_pairs_sql ?: [:]
        logMap += [default_hich_pairs_sql_opts: default_hich_pairs_sql_opts, hich_pairs_sql_opts: hich_pairs_sql_opts]
        def final_hich_pairs_sql_opts = buildCLIOpts(default_hich_pairs_sql_opts, hich_pairs_sql_opts, [:], null)
        def sql = hich_pairs_sql_opts.sql
        hich_pairs_sql_cmd = sql ? "hich pairs sql ${final_hich_pairs_sql_opts} '${sql}' /dev/stdin" : null
    }
    
    def default_pairtools_sort_opts = [
        "--output": output, 
        "--memory": "${memoryV}G", 
        "--nproc-in": cpus, 
        "--nproc-out": cpus,
        "--tmpdir": "."
    ]
    def pairtools_sort_opts = parse_to_pairs_opts?.pairtools_sort ?: [:]
    logMap += [default_pairtools_sort_opts: default_pairtools_sort_opts, pairtools_sort_opts: pairtools_sort_opts]
    def final_pairtools_sort_opts = buildCLIOpts(default_pairtools_sort_opts, pairtools_sort_opts, [:], null)

    def pairtools_sort_cmd = "pairtools sort ${final_pairtools_sort_opts}"   
    def cmd = [samtoolsViewCmd, pairtools_parse2_cmd, hich_pairs_sql_cmd, pairtools_sort_cmd].findAll{it}.join(" | ")
    logMap += [cmd: cmd]

    return [cmd, logMap, output]
}
