include { buildCLIOpts } from '../../util/cli.nf'

def buildCmd(id, pairs, fragmentIndex, select_pairs_opts, cpus) {
    def output = "${id}_select.pairs.gz"

    // 1. Separate user options into CLI args and filter criteria
    def userFilters = select_pairs_opts?.filters ?: [:]
    def userCliOpts = select_pairs_opts?.findAll { argName, _ -> argName != "filters" } ?: [:]

    // 2. Delegate complex logic to helper functions
    def filterExpression = _buildFilterExpression(userFilters, fragmentIndex as Boolean)
    def (preCmd, cliOpts) = _buildCliOpts(id, userCliOpts, userFilters.chroms, output, cpus)

    // 3. Assemble the final command string
    def cmdParts = [
        preCmd,
        "pairtools select",
        cliOpts,
        filterExpression,
        "'${pairs}'"
    ]
    def cmd = cmdParts.findAll { it }.join(" ")

    // 4. Create a structured log map for debugging
    def logMap = [
        task: "SELECT_PAIRS",
        cmd: cmd,
        output: [pairs: output],
        input: [id: id, pairs: pairs, select_pairs_opts: select_pairs_opts]
    ]

    return [cmd, logMap, output]
}

private String _buildFilterExpression(Map userFilters, boolean hasFragmentIndex) {
    // Establish default filter criteria
    def defaults = [
        discardSingleFrag: hasFragmentIndex,
        keepPairTypes: ["UU", "UR", "RU"]
    ]
    def filters = defaults + userFilters

    // Create a list of individual filter conditions
    def filterParts = [
        // Pair type filter (e.g., "pair_type in ['UU', 'UR', 'RU']")
        filters.keepPairTypes ? "pair_type in ['${filters.keepPairTypes.join("', '")}']" : null,

        // Cis/Trans filter (e.g., "chrom1 == chrom2")
        _getCisTransFilter(filters.onlyCis, filters.onlyTrans),

        // Strand-distance filter (e.g., "(strand1 + strand2 == '+-' and abs(pos2 - pos1) >= 1000)")
        _getStrandDistFilter(filters.minDistFR, filters.minDistRF, filters.minDistFF, filters.minDistRR),

        // Single fragment filter
        filters.discardSingleFrag ? "rfrag1 != rfrag2" : null,

        // User-provided custom filter string
        filters.custom
    ]

    // Join all valid conditions with "and"
    def finalFilter = filterParts.findAll { it }.collect { "(${it})" }.join(" and ") ?: "True"
    return "\"${finalFilter}\""
}

private List _buildCliOpts(String id, Map userCliOpts, List chroms, String output, int cpus) {
    // Modify --output-rest path to be sample-specific
    if ("--output-rest" in userCliOpts) {
        userCliOpts["--output-rest"] = "${id}${userCliOpts['--output-rest']}"
    }

    // Generate pre-command and CLI arg for --chrom-subset if needed
    def preCmd = ""
    def chromsFileOpt = ""
    if (chroms) {
        def bed = "__chroms__.bed"
        preCmd = "echo '${chroms.join('\\n')}' > '${bed}' &&"
        chromsFileOpt = "--chrom-subset '${bed}'"
    }

    def defaultOpts = ["--output": output, "--nproc-in": cpus, "--nproc-out": cpus]
    def remap = ["--output": "-o", "--type-cast": "-t", "--remove-columns": "-r"]
    def finalOpts = buildCLIOpts(defaultOpts, userCliOpts, remap, null)

    def allOpts = [chromsFileOpt, finalOpts].findAll { it }.join(" ")
    return [preCmd, allOpts]
}

private String _getCisTransFilter(Boolean onlyCis, Boolean onlyTrans) {
    if (onlyCis ^ onlyTrans) { // Use XOR to check if exactly one is true
        return onlyCis ? "chrom1 == chrom2" : "chrom1 != chrom2"
    }
    return null
}

private String _getStrandDistFilter(minDistFR, minDistRF, minDistFF, minDistRR) {
    def argMap = ["+-": minDistFR, "-+": minDistRF, "++": minDistFF, "--": minDistRR]
    def args = argMap.findAll { _, dist -> dist != null }.collect { strand, dist ->
        "(strand1 + strand2 == '${strand}' and abs(pos2 - pos1) >= ${dist})"
    }
    return args ? args.join(" or ") : null
}