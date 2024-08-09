include {parameterize} from './extraops.nf'

process HicrepCombos{
    publishDir "results/hicrep",
               mode: params.general.publish.mode
    container "bskubi/hich:latest"

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

    parameterize("hicrep",
                 samples,
                 params.comparisonSets,
                 ["mcool"],
                 ["mcool", "resolutions", "chroms", "exclude", "chrom_filter", "h", "dBPMax", "bDownSample"])
        | HicrepCombos

    emit:
    samples
}