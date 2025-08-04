include {validateMemory} from '../../util/memory.nf'
include {formatCLIArgs} from '../../util/cli.nf'
import java.lang.reflect.InvocationTargetException


def buildCmd(id, sambam, chromsizes, assembly, parseToPairs_opts, minMapq, memory, cpus) {
    def memoryV = validateMemory(memory, 2, 2)
    def output = "${id}.pairs.gz" 
    def samtoolsViewCmd = "samtools view -b '${sambam}'"

    def pairtoolsParse2_opts = ["--assembly": assembly, "--chroms-path": chromsizes, "--min-mapq": minMapq, "--flip": true]
    try{
        pairtoolsParse2_opts = formatCLIArgs(pairtoolsParse2_opts, parseToPairs_opts?.pairtools_parse2)
    } catch(InvocationTargetException e) {
        print(e.getTargetException())
    }
    def pairtoolsParse2Cmd = "pairtools parse2 ${pairtoolsParse2_opts}"
    
    def hichPairsSQLCmd = null
    if (parseToPairs_opts?.hichPairsSql?.sql) {
        def hichPairsSql_opts = ["--memory-limit": memoryV, "--threads": cpus]
        hichPairsSql_opts = formatCLIArgs(hichPairsSql_opts, parseToPairs_opts?.hich_pairs_sql)
        def sql = hichPairsSql_opts.sql
        hichPairsSQLCmd = sql ? "hich pairs sql ${hichPairsSql_opts} '${sql}' /dev/stdin" : null
    }
    
    def pairtoolsSort_opts = ["--output": output, "--memory": "${memoryV}G", "--nproc-in": cpus, "--nproc-out": cpus]
    pairtoolsSort_opts = formatCLIArgs(pairtoolsSort_opts, parseToPairs_opts?.pairtools_sort)

    def pairtoolsSortCmd = "pairtools sort ${pairtoolsSort_opts}"   
    def cmd = [samtoolsViewCmd, pairtoolsParse2Cmd, hichPairsSQLCmd, pairtoolsSortCmd].findAll{it}.join(" | ")

    def logMap = [
        task: "PARSE",
        input: [
            id: id,
            sambam: sambam,
            chromsizes: chromsizes,
            assembly: assembly,
            parseToPairs_opts: parseToPairs_opts,
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
