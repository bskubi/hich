include {buildCLIOpts} from '../../util/cli.nf'

def buildCmd(id, sambam, chromsizes, assembly, parse_to_pairs_opts, minMapq, memory, cpus) {
    def output = "${id}.pairs.gz"

    // 1. Call helper functions to build each part of the pipeline
    def samtoolsCmd = _buildSamtoolsCmd(sambam)
    def parseResult = _buildPairtoolsParseCmd(assembly, chromsizes, minMapq, parse_to_pairs_opts)
    def sqlResult   = _buildHichSqlCmd(memory, cpus, parse_to_pairs_opts)
    def sortResult  = _buildPairtoolsSortCmd(output, memory, cpus, parse_to_pairs_opts)

    // 2. Assemble the final command from the parts, skipping any nulls (like the optional sql step)
    def cmdParts = [
        samtoolsCmd,
        parseResult.cmd,
        sqlResult?.cmd, // Safe navigation for the optional command
        sortResult.cmd
    ]
    def cmd = cmdParts.findAll { it }.join(" | ")

    // 3. Assemble log map
    def logMap = [
        task: "PARSE_TO_PAIRS",
        cmd: cmd,
        output: [pairs: output],
        input: [id: id, sambam: sambam, chromsizes: chromsizes, assembly: assembly],
        opts_used: [
            pairtools_parse2: [defaults: parseResult.defaults, user: parseResult.user],
            hich_pairs_sql:   [defaults: sqlResult?.defaults, user: sqlResult?.user],
            pairtools_sort:   [defaults: sortResult.defaults, user: sortResult.user]
        ]
    ]

    return [cmd, logMap, output]
}

private String _buildSamtoolsCmd(String sambam) {
    return "samtools view -b '${sambam}'"
}

private Map _buildPairtoolsParseCmd(String assembly, String chromsizes, int minMapq, Map opts) {
    def defaultOpts = [
        "--assembly": assembly,
        "--chroms-path": chromsizes,
        "--min-mapq": minMapq,
        "--flip": true
    ]
    def userOpts = opts?.pairtools_parse2 ?: [:]
    def remap = ["--chroms-path": "-c"]
    def finalOpts = buildCLIOpts(defaultOpts, userOpts, remap, null)

    return [
        cmd: "pairtools parse2 ${finalOpts}",
        defaults: defaultOpts,
        user: userOpts
    ]
}

private Map _buildHichSqlCmd(int memory, int cpus, Map opts) {
    def hichSqlOpts = opts?.hich_pairs_sql
    if (!hichSqlOpts?.sql) {
        return null // This command is optional
    }

    def defaultOpts = ["--memory-limit": memory, "--threads": cpus]
    def finalOpts = buildCLIOpts(defaultOpts, hichSqlOpts, [:], null)

    return [
        cmd: "hich pairs sql ${finalOpts} '${hichSqlOpts.sql}' /dev/stdin",
        defaults: defaultOpts,
        user: hichSqlOpts
    ]
}

private Map _buildPairtoolsSortCmd(String output, int memory, int cpus, Map opts) {
    def defaultOpts = [
        "--output": output,
        "--memory": "${memory}G",
        "--nproc-in": cpus,
        "--nproc-out": cpus,
        "--tmpdir": "."
    ]
    def userOpts = opts?.pairtools_sort ?: [:]
    def finalOpts = buildCLIOpts(defaultOpts, userOpts, [:], null)

    return [
        cmd: "pairtools sort ${finalOpts}",
        defaults: defaultOpts,
        user: userOpts
    ]
}

