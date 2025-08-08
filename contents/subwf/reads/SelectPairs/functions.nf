include {buildCLIOpts} from '../../util/cli.nf'

def dq(it) {
    "\"${it}\""
}

def sq(it) {
    "'${it}'"
}

def sb(it) {
    "[${it}]"
}

def pythonList(it) {
    it = it.collect{sq(it)}
    it = it.join(", ")
    it = it ? sb(it) : null
}

def formatKeepPairTypes(pairTypes) {
    def keep = pythonList(pairTypes)
    def arg = keep ? "pair_type in ${keep}" : null
    return arg
}

def formatCisTrans(onlyCis, onlyTrans) {
    onlyCis = onlyCis ?: false
    onlyTrans = onlyTrans ?: false
    def arg = null
    if (onlyCis ^ onlyTrans) {
        arg = onlyCis ? "chrom1 == chrom2" : "chrom1 != chrom2"
    }
    return arg
}

def formatStrandDistFilters(minDistFR, minDistRF, minDistFF, minDistRR) {
    def argMap = ["+-": minDistFR, "-+": minDistRF, "++": minDistFF, "--": minDistRR]
    def args = []
    argMap.each {
        strand, dist ->
        if (dist != null) {
            strand = sq(strand)
            args = args + ["(strand1 + strand2 == ${strand} and abs(pos2 - pos1) >= ${dist})"]
        }
    }
    args = args.join(" or ")
    args = args ?: null
    return args
}

def formatDiscardSingleFrag(discardSingleFrag) {
    return discardSingleFrag ? "rfrag1 != rfrag2" : null
}

def formatFilters(filters) {
    filters = filters.findAll{it}
    filters = filters.collect {"(${it})"}
    filters = filters.join(" and ")
    filters = filters ?: "True"
    filters = dq(filters)
    return filters
}

def updateOutputPaths(id, pairtoolsSelect_opts) {
    if ("--output-rest" in pairtoolsSelect_opts) {
        def rest = pairtoolsSelect_opts["--output-rest"]
        rest = "${id}_${rest}"
        pairtoolsSelect_opts += ["--output-rest": rest]
    }
    return pairtoolsSelect_opts
}

def formatWriteChroms(chroms) {
    chroms = chroms ?: []
    chroms = chroms.join('\n')
    chroms = chroms ? sq(chroms) : ""
    def bed = sq("__chroms__.bed")
    chroms = chroms ? "echo ${chroms} > ${bed} &&" : ""
    def chromsFile = chroms ? "--chrom-subset ${bed}" : ""
    return [chroms, chromsFile]
}

def buildPairtoolsSelectFilters(default_filters, user_filters) {
    def filters = (default_filters ?: [:]) + (user_filters ?: [:])
    def pairTypes = formatKeepPairTypes(filters.keepPairTypes)
    def cisTrans = formatCisTrans(filters.onlyCis, filters.onlyTrans)
    def strandDist = formatStrandDistFilters(filters.minDistFR, filters.minDistRF, filters.minDistFF, filters.minDistRR)
    def discardSingleFrag = formatDiscardSingleFrag(filters.discardSingleFrag)
    def custom = filters.custom
    filters = formatFilters([pairTypes, cisTrans, strandDist, discardSingleFrag, custom])

    return filters
}

def buildCmd(id, pairs, fragmentIndex, select_pairs_opts, cpus) {
    def output = "${id}_select.pairs.gz"
    def logMap = [
        task: "SELECT_PAIRS",
        input: [
            id: id,
            pairs: pairs,
            select_pairs_opts: select_pairs_opts
        ],
        output: [
            pairs: output
        ]
    ]
    def pairtools_select_filters = select_pairs_opts?.filters ?: [:]
    def pairtools_select_opts = select_pairs_opts?.findAll{argName, argVal -> argName != "filters"}
    pairtools_select_opts = updateOutputPaths(id, pairtools_select_opts)

    def default_pairtools_select_filters = [
        discardSingleFrag: fragmentIndex as Boolean,
        keepPairTypes: ["UU", "UR", "RU"]
    ]
    def filters = buildPairtoolsSelectFilters(default_pairtools_select_filters, pairtools_select_filters)
    def (write_chroms, chroms_file) = formatWriteChroms(pairtools_select_filters.chroms)
    
    def default_pairtools_select_opts = ["--output": output, "--nproc-in": cpus, "--nproc-out": cpus]

    def remap = [
        "--output": "-o",
        "--type-cast": "-t",
        "--remove-columns": "-r"
    ]
    def final_pairtools_select_opts = buildCLIOpts(default_pairtools_select_opts, pairtools_select_opts, remap, null)
    def pairsInput = sq(pairs)

    def cmd = [
        write_chroms,
        "pairtools select",
        chroms_file,
        final_pairtools_select_opts,
        filters,
        pairsInput
    ]
    cmd = cmd.findAll{it}
    cmd = cmd.join(" ")
    logMap.cmd = cmd

    return [cmd, logMap, output]
}
