process MustacheDiffloops{
    publishDir "results/loops"
    container "mustache"

    input:
    tuple val(prefix), val(id1), path(mx1), val(id2), path(mx2), val(mustache_params)

    output:
    tuple val(id1), val(id2), path("${prefix}.loop1"), path("${prefix}.loop2"), path("${prefix}.diffloop1"), path("${prefix}.diffloop2")

    shell:
    cmd = ["python /mustache/mustache/diff_mustache.py",
           "-f1 ${mx1} -f2 ${mx2}",
           "-o ${prefix}"] + mustache_params
    cmd = cmd.join(" ")
    cmd

    stub:
    "touch ${prefix}.loop1 ${prefix}.loop2 ${prefix}.diffloop1 ${prefix}.diffloop2"
}

workflow CallLoops {
    take:
    samples

    main:
    // For feature calling, we have to segregate by merge (techrep/biorep/condition)
    // and by whether or not it is downsampled

    samples
        | filter {
            sample ->
            (sample.is_condition && "is_condition" in sample.loops.call_on) ||
            (!sample.is_condition && sample.is_biorep && "is_biorep" in sample.loops.call_on) ||
            (!sample.is_condition && !sample.is_biorep && sample.is_techrep && "is_techrep" in sample.loops.call_on)
        }
        | map{sample -> sample.subMap("id", "mcool", "hic", "is_techrep", "is_biorep", "is_condition", "loops")}
        | map{[[it.subMap("is_techrep", "is_biorep", "is_condition")], it]}
        | groupTuple
        | map{it[1]}
        | map {
            comparison_set ->
            comparisons = []
            comparison_set.eachWithIndex {
                sample1, idx1->
                sub_list = idx1 + 1 < comparison_set.size() ? comparison_set[(idx1+1)..-1] : []
                sub_list.each {
                    sample2 ->
                    comparisons += [[sample1, sample2]]
                }
            }
            comparisons
        }
        | flatMap
        | map{["${it[0].id}_${it[1].id}", it[0].id, it[0].mcool, it[1].id, it[1].mcool, it[0].loops.mustache_params]}
        | MustacheDiffloops

    emit:
    samples
}