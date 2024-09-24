include {parameterize} from './extraops.nf'

process HicrepCombos{
    publishDir "results/hicrep",
               mode: params.general.publish.mode
    //container "bskubi/hich:latest"

    input:
    tuple path(mcools), val(resolutions), val(chroms), val(exclude), val(chromFilter), val(h), val(dBPMax), val(bDownSample)

    output:
    path("hicrep.tsv")

    shell:

    cmd = ["hich hicrep",
           resolutions ? "--resolutions ${resolutions.join(",")}" : "",
           chroms ? "--chroms ${chroms.join(",")}" : "",
           exclude ? "--exclude ${exclude.join(",")}" : "",
           chromFilter ? "--chrom-filter '${chromFilter}'" : "",
           "--h ${h.join(",")}",
           "--d-bp-max ${dBPMax.join(",")}",
           "--b-downsample ${bDownSample.join(",")}",
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
    
    parameterize("hicrep",
                 samples,
                 params.comparisonSets,
                 ["mcool"],
                 ["mcool", "resolutions", "chroms", "exclude", "chromFilter", "h", "dBPMax", "bDownSample"])
        | HicrepCombos

    emit:
    samples
}