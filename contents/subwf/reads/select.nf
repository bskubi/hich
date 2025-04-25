include {QCReads} from './qcHicReads.nf'
include {emptyOnLastStep; pack; skip} from '../extraops.nf'
include {withLog; stubLog} from '../util/logs.nf'

process PairtoolsSelect {
    publishDir params.general.publish.select ? params.general.publish.select : "results",
               saveAs: {params.general.publish.select ? it : null},
               mode: params.general.publish.mode
    label 'doJobArray'
    label 'pairs'
    cpus 8

    input:
    tuple val(id), path(pairs), val(pairtoolsSelectParams), val(condition)

    output:
    tuple val(id), path("${id}_select.pairs.gz")

    shell:
    def quote = { List<String> items ->
    result = items.collect { "'$it'" }.join(", ")
    "[${result}]"}

    pair_types = "pair_type in ${quote(condition.keepPairTypes)}"
    
    cis_trans = condition.keepCis ^ condition.keepTrans
                    ? (condition.keepCis
                            ? "(chrom1 == chrom2)"
                            : "(chrom1 != chrom2)")
                    : null

    min_distances = [:]
    min_distances += condition.minDistFR != null ? ["+-":condition.minDistFR] : [:]
    min_distances += condition.minDistRF != null ? ["-+":condition.minDistRF] : [:]
    min_distances += condition.minDistFF != null ? ["++":condition.minDistFF] : [:]
    min_distances += condition.minDistRR != null ? ["--":condition.minDistRR] : [:]
    strand_dist = min_distances.collect {
        strand, dist ->
        s1 = strand[0]
        s2 = strand[1]
        "(strand1 + strand2 == '${s1}${s2}' and abs(pos2 - pos1) >= ${dist})"
    }.join(" or ")

    strand_dist = strand_dist ?: null
    
    frags = condition.discardSingleFrag ? "rfrag1 != rfrag2" : null
    
    filters = [pair_types, cis_trans, strand_dist, frags]
    
    filters.removeAll([null])
    filters = filters.collect {"(${it})"}
    filters = filters.join(" and ")

    pairtoolsSelectParams = pairtoolsSelectParams ? pairtoolsSelectParams.collect {
        item ->
        ["--output-rest": "--output-rest ${id}_select.rest.pairs.gz"].get(item, item)
    }.join(" ") : ""

    write_chroms = condition.chroms ? "echo '${condition.chroms.join('\n')}' > __chroms__.bed &&" : ""
    chroms_file = condition.chroms ? "--chrom-subset __chroms__.bed" : ""

    cmd = """${write_chroms} pairtools select --output '${id}_select.pairs.gz' ${chroms_file} ${pairtoolsSelectParams} "${filters}"  --nproc-in ${task.cpus} --nproc-out ${task.cpus} '${pairs}'"""
    
    logMap = [
        task: "PairtoolsSelect",
        input: [
            id: id,
            pairs: pairs,
            pairtoolsSelectParams: pairtoolsSelectParams,
            condition: condition
        ],
        output: [
            pairs: "${id}_select.pairs.gz"
        ]
    ]

    withLog(cmd, logMap)


    stub:
    stub = "touch '${id}_select.pairs.gz'"
    def quote = { List<String> items ->
    result = items.collect { "'$it'" }.join(", ")
    "[${result}]"}

    pair_types = "pair_type in ${quote(condition.keepPairTypes)}"
    
    cis_trans = condition.keepCis ^ condition.keepTrans
                    ? (condition.keepCis
                            ? "(chrom1 == chrom2)"
                            : "(chrom1 != chrom2)")
                    : null

    min_distances = [:]
    min_distances += condition.minDistFR != null ? ["+-":condition.minDistFR] : [:]
    min_distances += condition.minDistRF != null ? ["-+":condition.minDistRF] : [:]
    min_distances += condition.minDistFF != null ? ["++":condition.minDistFF] : [:]
    min_distances += condition.minDistRR != null ? ["--":condition.minDistRR] : [:]
    strand_dist = min_distances.collect {
        strand, dist ->
        s1 = strand[0]
        s2 = strand[1]
        "(strand1 + strand2 == '${s1}${s2}' and abs(pos2 - pos1) >= ${dist})"
    }.join(" or ")

    strand_dist = strand_dist ?: null
    
    frags = condition.discardSingleFrag ? "rfrag1 != rfrag2" : null
    
    filters = [pair_types, cis_trans, strand_dist, frags]
    
    filters.removeAll([null])
    filters = filters.collect {"(${it})"}
    filters = filters.join(" and ")

    pairtoolsSelectParams = pairtoolsSelectParams ? pairtoolsSelectParams.collect {
        item ->
        ["--output-rest": "--output-rest ${id}_select.rest.pairs.gz"].get(item, item)
    }.join(" ") : ""

    write_chroms = condition.chroms ? "echo '${condition.chroms.join('\n')}' > __chroms__.bed &&" : ""
    chroms_file = condition.chroms ? "--chrom-subset __chroms__.bed" : ""

    cmd = """${write_chroms} pairtools select --output '${id}_select.pairs.gz' ${chroms_file} ${pairtoolsSelectParams} "${filters}"  --nproc-in ${task.cpus} --nproc-out ${task.cpus} '${pairs}'"""
    
    logMap = [
        task: "PairtoolsSelect",
        input: [
            id: id,
            pairs: pairs,
            pairtoolsSelectParams: pairtoolsSelectParams,
            condition: condition
        ],
        output: [
            pairs: "${id}_select.pairs.gz"
        ]
    ]

    stubLog(stub, cmd, logMap)
}

workflow Select {
    take:
    samples
    
    main:

    samples
        | filter{!skip("select") && (it.pairtoolsSelectParams || it.pairtoolsSelectFilters) && (it.pairs || it.latestPairs)}
        | map{tuple(it.id, it.latestPairs, it.pairtoolsSelectParams, it.pairtoolsSelectFilters)}
        | PairtoolsSelect
        | map{[id:it[0], selectPairs:it[1], latest:it[1], latestPairs:it[1]]}
        | set{result}
    pack(samples, result) | set{samples}

    
    if ("select" in params.general.get("qcAfter")) {
        QCReads(samples, "select")
    }

    samples = emptyOnLastStep("select", samples)




    emit:
    samples
}