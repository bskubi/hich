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

def formatPairtoolsSelectParams(id, pairtoolsSelectParams) {
    def rest = sq("${id}_select.rest.pairs.gz")
    pairtoolsSelectParams = pairtoolsSelectParams ? pairtoolsSelectParams.collect {
            item ->
            
            ["--output-rest": "--output-rest ${rest}"].get(item, item)
        }.join(" ") : ""
    return pairtoolsSelectParams
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

def buildCmd(id, pairs, pairtoolsSelectParams, pairtoolsSelectFilters, cpus) {
    def pairTypes = formatKeepPairTypes(pairtoolsSelectFilters.keepPairTypes)
    def cisTrans = formatCisTrans(pairtoolsSelectFilters.onlyCis, pairtoolsSelectFilters.onlyTrans)
    def strandDist = formatStrandDistFilters(pairtoolsSelectFilters.minDistFR, pairtoolsSelectFilters.minDistRF, pairtoolsSelectFilters.minDistFF, pairtoolsSelectFilters.minDistRR)
    def discardSingleFrag = formatDiscardSingleFrag(pairtoolsSelectFilters.discardSingleFrag)
    def filters = formatFilters([pairTypes, cisTrans, strandDist, discardSingleFrag])
    pairtoolsSelectParams = formatPairtoolsSelectParams(id, pairtoolsSelectParams)
    def (writeChroms, chromsFile) = formatWriteChroms(pairtoolsSelectFilters.chroms)
    def output = "${id}_select.pairs.gz"
    def pairsOutput = "--output ${sq(output)}"
    def nprocIn = "--nproc-in ${cpus}"
    def nprocOut = "--nproc-out ${cpus}"
    def pairsInput = sq(pairs)



    def cmd = [
        writeChroms,
        "pairtools select",
        pairsOutput,
        chromsFile,
        pairtoolsSelectParams,
        nprocIn,
        nprocOut,
        filters,
        pairsInput
    ]
    cmd = cmd.findAll{it}
    cmd = cmd.join(" ")
    
    logMap = [
        task: "PairtoolsSelect",
        input: [
            id: id,
            pairs: pairs,
            pairtoolsSelectParams: pairtoolsSelectParams,
            pairtoolsSelectFilters: pairtoolsSelectFilters
        ],
        output: [
            pairs: output
        ]
    ]

    return [cmd, logMap, output]
}
