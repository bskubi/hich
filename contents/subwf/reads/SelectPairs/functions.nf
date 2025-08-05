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

def buildFilters(filters) {
    filters = filters ?: [:]
    def pairTypes = formatKeepPairTypes(filters.keepPairTypes)
    def cisTrans = formatCisTrans(filters.onlyCis, filters.onlyTrans)
    def strandDist = formatStrandDistFilters(filters.minDistFR, filters.minDistRF, filters.minDistFF, filters.minDistRR)
    def discardSingleFrag = formatDiscardSingleFrag(filters.discardSingleFrag)
    filters = formatFilters([pairTypes, cisTrans, strandDist, discardSingleFrag])

    return filters
}

def buildCmd(id, pairs, selectPairs_opts, cpus) {
    def pairtoolsSelectFilters = selectPairs_opts?.filters ?: [:]
    def pairtoolsSelect_opts = selectPairs_opts?.findAll{argName, argVal -> argName != "filters"}
    pairtoolsSelect_opts = updateOutputPaths(id, pairtoolsSelect_opts)
    def filters = buildFilters(pairtoolsSelectFilters)
    def (writeChroms, chromsFile) = formatWriteChroms(pairtoolsSelectFilters.chroms)

    def output = "${id}_select.pairs.gz"
    def defaultPairtoolsSelect_opts = ["--output": output, "--nproc-in": cpus, "--nproc-out": cpus]
    pairtoolsSelect_opts = buildCLIOpts(defaultPairtoolsSelect_opts, pairtoolsSelect_opts)
    def pairsInput = sq(pairs)

    def cmd = [
        writeChroms,
        "pairtools select",
        chromsFile,
        pairtoolsSelect_opts,
        filters,
        pairsInput
    ]
    cmd = cmd.findAll{it}
    cmd = cmd.join(" ")

    logMap = [
        task: "SELECT_PAIRS",
        input: [
            id: id,
            pairs: pairs,
            selectPairs_opts: selectPairs_opts
        ],
        output: [
            pairs: output
        ]
    ]

    return [cmd, logMap, output]
}
