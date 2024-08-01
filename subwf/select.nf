include {QCReads} from './qcHicReads.nf'
include {transpack} from './extraops.nf'

process PairtoolsSelect {
    publishDir params.general.publish.select ? params.general.publish.select : "results",
               saveAs: {params.general.publish.select ? it : null},
               mode: params.general.publish.mode
    container "bskubi/hich:latest"

    input:
    tuple val(id), path(pairs), val(select_params), val(condition)

    output:
    tuple val(id), path("${id}_select.pairs.gz")

    shell:
    def quote = { List<String> items ->
    result = items.collect { "'$it'" }.join(", ")
    "[${result}]"}

    pair_types = "pair_type in ${quote(condition.keep_pair_types)}"
    
    cis_trans = condition.keep_cis ^ condition.keep_trans
                    ? (condition.keep_cis
                            ? "(chrom1 == chrom2)"
                            : "(chrom1 != chrom2)")
                    : null

    min_distances = [:]
    min_distances += condition.min_dist_fr != null ? ["+-":condition.min_dist_fr] : [:]
    min_distances += condition.min_dist_rf != null ? ["-+":condition.min_dist_rf] : [:]
    min_distances += condition.min_dist_ff != null ? ["++":condition.min_dist_ff] : [:]
    min_distances += condition.min_dist_rr != null ? ["--":condition.min_dist_rr] : [:]
    strand_dist = min_distances.collect {
        strand, dist ->
        s1 = strand[0]
        s2 = strand[1]
        "(strand1 + strand2 == '${s1}${s2}' and abs(pos2 - pos1) >= ${dist})"
    }.join(" or ")

    strand_dist = strand_dist ?: null
    
    frags = condition.discard_same_frag ? "rfrag1 != rfrag2" : null
    
    filters = [pair_types, cis_trans, strand_dist, frags]
    
    filters.removeAll([null])
    filters = filters.collect {"(${it})"}
    filters = filters.join(" and ")

    select_params = select_params ? select_params.collect {
        item ->
        ["--output-rest": "--output-rest ${id}_select.rest.pairs.gz"].get(item, item)
    }.join(" ") : ""

    write_chroms = condition.chroms ? "echo '${condition.chroms.join('\n')}' > __chroms__.bed &&" : ""
    chroms_file = condition.chroms ? "--chrom-subset __chroms__.bed" : ""

    cmd = """${write_chroms} pairtools select --output ${id}_select.pairs.gz ${chroms_file} ${select_params} "${filters}" ${pairs}"""
    
    cmd

    stub:
    "touch ${id}_select.pairs.gz"
}

workflow Select {
    take:
        samples
    
    main:        
        transpack (
            PairtoolsSelect,
            samples,
            ["id", "latest", "select_params", "select_condition"],
            ["id", "select_pairs"],
            ["latest":"select_pairs"],
            "id",
            ["nullOk":"select_params"]
        ) | set{samples}
        
        if ("Select" in params.general.get("qc_after")) {
            QCReads(samples, "Select")
        }

        if (params.general.get("last_step") == "select") {
            channel.empty() | set{samples}
        }




    emit:
        samples
}