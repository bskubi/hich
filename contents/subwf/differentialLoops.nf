include {label; isTechrep; isBiorep; isCondition; parameterize; combinations} from './extraops.nf'

process MustacheDiffloops{
    publishDir "results/loops",
               mode: params.general.publish.mode
    container "bskubi/mustache:latest"

    input:
    tuple val(prefix), path(mx1), path(mx2), val(mustache_params)

    output:
    tuple path("${prefix}.loop1"), path("${prefix}.loop2"), path("${prefix}.diffloop1"), path("${prefix}.diffloop2")

    shell:
    cmd = ["python /mustache/mustache/diff_mustache.py",
           "-f1 ${mx1} -f2 ${mx2}",
           "-o ${prefix}"] + mustache_params
    cmd = cmd.join(" ")
    cmd

    stub:
    "touch ${prefix}.loop1 ${prefix}.loop2 ${prefix}.diffloop1 ${prefix}.diffloop2"
}

process MustacheLoops{
    publishDir "results/loops",
               mode: params.general.publish.mode
    container "bskubi/mustache:latest"

    input:
    tuple val(prefix), path(mx), val(mustache_params)

    output:
    path("${prefix}.loop")

    shell:
    cmd = ["python -m mustache -f ${mx} -o ${prefix}.loop"] + mustache_params
    cmd = cmd.join(" ")
    cmd

    stub:
    "touch ${prefix}.loop1 ${prefix}.loop2 ${prefix}.diffloop1 ${prefix}.diffloop2"
}

workflow DifferentialLoops {
    take:
    samples

    main:

    /*
        We need to alter this so that the user can control:
            Which aggregation profiles to run it on
            Which aggregation levels to compare
            Which samples to exclude from comparisons
        
        This might look like:
            
    */

    samples
        | branch {
            sample ->

            techreps: isTechrep(sample)
            bioreps: isBiorep(sample)
            conditions: isCondition(sample)
        } | set{branched}
    
    mustacheDiffloopsInput = channel.empty()

    [branched.techreps, branched.bioreps, branched.conditions].each {
        sampleTypeGroupChannel ->
        p = parameterize("mustacheLoops",
                         sampleTypeGroupChannel,
                         params.comparisonSets,
                         ["id", "latestMatrix"],
                         ["id", "latestMatrix", "mustacheParams"])

        combinations(p,
                     ["id", "latestMatrix"],
                     ["id1", "mx1", "id2", "mx2"],
                     ["mustacheParams"])
            | map{tuple("${it.id1}_${it.id2}".toString(), it.mx1, it.mx2, it.mustacheParams)}
            | set{c}
        
        mustacheDiffloopsInput = mustacheDiffloopsInput.concat(c)
    }

    mustacheDiffloopsInput | MustacheDiffloops

    emit:
    samples
}