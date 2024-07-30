process HicrepCombos{
    publishDir "results/hicrep", mode:"copy"
    container "hich"

    input:
    tuple path(mcools), val(resolutions), val(chroms), val(exclude), val(chrom_filter), val(h), val(dBPMax), val(bDownSample)

    output:
    path("hicrep.tsv")

    shell:

    cmd = ["hich hicrep",
           resolutions ? "--resolutions ${resolutions.join(",")}" : "",
           chroms ? "--chroms ${chroms.join(",")}" : "",
           exclude ? "--exclude ${exclude.join(",")}" : "",
           chrom_filter ? "--chrom_filter '${chrom_filter}'" : "",
           "--h ${h.join(",")}",
           "--d_bp_max ${dBPMax.join(",")}",
           "--b_downsample ${bDownSample.join(",")}",
           "--output hicrep.tsv",
           "${mcools.join(" ")}"].join(" ")
    cmd

    stub:
    "touch hicrep.tsv"
}

workflow Hicrep {
    take:
    samples

    main:

    samples
        | filter {
            sample ->
            (sample.is_condition && "is_condition" in sample.hicrep.call_on) ||
            (!sample.is_condition && sample.is_biorep && "is_biorep" in sample.hicrep.call_on) ||
            (!sample.is_condition && !sample.is_biorep && sample.is_techrep && "is_techrep" in sample.hicrep.call_on)
        }
        | filter{sample -> sample.containsKey("mcool")}
        | map{sample -> sample.get("mcool")}
        | collect
        | map{mcool -> tuple(mcool,
                             params.defaults.hicrep.resolutions,
                             params.defaults.hicrep.chroms,
                             params.defaults.hicrep.exclude,
                             params.defaults.hicrep.chrom_filter,
                             params.defaults.hicrep.h,
                             params.defaults.hicrep.dBPMax,
                             params.defaults.hicrep.bDownSample)}
        | HicrepCombos

    emit:
    samples
}