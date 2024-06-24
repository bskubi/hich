process align_hic {
    tag       "${sample_id}"
    container "bskubi/bwa-mem2:latest"
    maxRetries 4
    memory {20.GB + 10.GB * task.attempt}

    //Alignment speed scales ~linearly with processes, so does not benefit
    //from running alignment on multiple files at once. Each instance of
    //bwa-mem2 uses about 15 gigs of memory. For these two reasons we 
    //tell nextflow to run one alignment process at a time with maxForks 1.
    maxForks 1

    input:
        tuple val(sample_id), path(reads, arity: 2), path(index_dir), val(index_prefix)

    output:
        tuple val(sample_id), path("${sample_id}.bam")

    shell:
        [
            "bwa-mem2 mem -SP5M -t 36",
                "!{index_dir}/!{index_prefix}",
                "!{reads[0]} !{reads[1]}",
            "| samtools view -b -o !{sample_id}.bam"
        ].join(" ")
    
    stub:
        "touch ${sample_id}.bam"
}

def format_align_hic_input = {tuple(it.sample_id, tuple(it.fastq1, it.fastq2), file(it.index_dir), it.index_prefix)}
def format_align_hic_output = {[sample_id: it[0], bam: it[1]]}

process parse_to_contacts {
    publishDir "results/parse"
    tag       "${sample_id}"
    container "bskubi/pairtools:1.0.4"
    memory 4.Gb

    input:
        tuple val(sample_id), path(bam), path(chromsizes), val(flip), val(assembly), val(min_mapq), val(max_inter_align_gap), val(max_insert_size), val(dedup_max_mismatch)

    output:
        tuple val(sample_id), path("${sample_id}.parse2.pairs.gz"), path("${sample_id}.parse2.stats.txt")

    shell:
        start_cmd =[
            "pairtools parse2",
                "--chroms-path !{chromsizes}",
                "--output-stats !{sample_id}.parse2.stats.txt",
                "--output !{sample_id}.parse2.pairs.gz"
        ]
        optional_params = [
            assembly ? "--assembly !{assembly}" : "",
            flip != null ? "--flip" : "",
            min_mapq != null ? "--min-mapq !{min_mapq}" : "",
            max_inter_align_gap != null ? "--max_inter_align_map !{max_inter_align_gap}" : "",
            max_insert_size != null ? "--max_insert_size !{max_insert_size}" : "",
            dedup_max_mismatch != null ? "--dedup_max_mismatch !{dedup_max_mismatch}" : ""
        ].findAll{it}
        final_params = [
            "!{bam}"
        ]

        cmd = (start_cmd + optional_params + final_params).join(" ")
        cmd

    stub:
        "touch ${sample_id}.parse2.pairs.gz && touch ${sample_id}.parse2.stats.txt"
}

process restrict {
    tag       "${sample_id}"
    memory 5.Gb
    conda 'python-intervals lz4 hich-restrict=0.1.6'

    input:
    tuple( val(sample_id), val(frags_format), path(frags), path(pairs))

    output:
    tuple( val(sample_id), path("${sample_id}.fragtag.pairs.gz"))

    shell:
        cmd = [
            "hich-restrict --help > help_output.txt;",
            "hich-restrict",
            frags_format ? "--frags-file-format !{frags_format}" : "",
            "--output-pairs !{sample_id}.fragtag.pairs.gz",
            "!{frags}",
            "!{pairs}",
	    "&& echo done"
        ].findAll{it}.join(" ")

        cmd
    
    stub:
        "touch ${sample_id}.fragtag.pairs.gz"
}

def assemble_filter_contacts_proc_cmd(sample_id,
                                      pairs,
                                      select_map_type,
                                      filter_below_FF,
                                      filter_below_RR,
                                      filter_below_FR,
                                      filter_below_RF,
                                      filter_same_frag) 
{
    distance = "abs(pos2 - pos1)"
    strands = "(strand1, strand2)"

    if (!select_map_type) {
        select_map_type = ["UU", "UR", "RU"]
    }

    map_type_filter = "(" + select_map_type.collect{"""(pair_type=="${it}")"""}.join(" or ") + ")"

    sdf_in_use = [filter_below_FF, filter_below_RR, filter_below_FR, filter_below_RF].any{it != null}

    strand_distance_filter = sdf_in_use
        ? "(" + [filter_below_FF != null ? """(${strands} == ("+", "+") and ${distance} >= ${filter_below_FF})""" : """(${strands} == ("+", "+"))""",
                 filter_below_RR != null ? """(${strands} == ("-", "-") and ${distance} >= ${filter_below_RR})""" : """(${strands} == ("-", "-"))""",
                 filter_below_FR != null ? """(${strands} == ("+", "-") and ${distance} >= ${filter_below_FR})""" : """(${strands} == ("+", "-"))""",
                 filter_below_RF != null ? """(${strands} == ("-", "+") and ${distance} >= ${filter_below_RF})""" : """(${strands} == ("-", "+"))""",
            ].findAll{it}.join(" or ") + ")"
        : ""
    
    same_frag_filter = filter_same_frag ? """(COLS[-6] != COLS[-3])""" : ""

    all_filters = [
        map_type_filter,
        strand_distance_filter,
        same_frag_filter
    ].findAll{it}.join(" and ")

    filter_cmd = all_filters ? "${all_filters}" : "True"

    cmd = [
        "pairtools select",
        "'${filter_cmd}'",
        "--output-rest ${sample_id}.filtered_out.pairs.gz",
        "--output ${sample_id}.select.pairs.gz",
        "${pairs}"

    ].join(" ")
    return cmd
}

process filter_contacts {
    publishDir "results/filter"
    tag       "${sample_id}"
    container "bskubi/pairtools:1.0.4"
    memory 4.Gb

    input:
        tuple(val(sample_id),
              path(pairs),
              val(select_map_type),
              val(filter_below_FF),
              val(filter_below_RR),
              val(filter_below_FR),
              val(filter_below_RF),
              val(filter_same_frag))

    output:
        tuple val(sample_id), path("${sample_id}.select.pairs.gz")

    shell:
        distance = "abs(pos2 - pos1)"
        strands = "(strand1, strand2)"

        if (!select_map_type) {
            select_map_type = ["UU", "UR", "RU"]
        }

        map_type_filter = "(" + select_map_type.collect{"""(pair_type=="${it}")"""}.join(" or ") + ")"

        sdf_in_use = [filter_below_FF, filter_below_RR, filter_below_FR, filter_below_RF].any{it != null}

        strand_distance_filter = sdf_in_use
            ? "(" + [filter_below_FF != null ? """(${strands} == ("+", "+") and ${distance} >= ${filter_below_FF})""" : """(${strands} == ("+", "+"))""",
                    filter_below_RR != null ? """(${strands} == ("-", "-") and ${distance} >= ${filter_below_RR})""" : """(${strands} == ("-", "-"))""",
                    filter_below_FR != null ? """(${strands} == ("+", "-") and ${distance} >= ${filter_below_FR})""" : """(${strands} == ("+", "-"))""",
                    filter_below_RF != null ? """(${strands} == ("-", "+") and ${distance} >= ${filter_below_RF})""" : """(${strands} == ("-", "+"))""",
                ].findAll{it}.join(" or ") + ")"
            : ""
        
        same_frag_filter = filter_same_frag ? """(COLS[-6] != COLS[-3])""" : ""

        all_filters = [
            map_type_filter,
            strand_distance_filter,
            same_frag_filter
        ].findAll{it}.join(" and ")

        filter_cmd = all_filters ? "${all_filters}" : "True"

        cmd = [
            "pairtools select",
            "'${filter_cmd}'",
            "--output-rest ${sample_id}.filtered_out.pairs.gz",
            "--output ${sample_id}.select.pairs.gz",
            "${pairs}"

        ].join(" ")
        println("\n${sample_id} ${pairs} ${cmd}\n")
        cmd

    stub:
        // println(assemble_filter_contacts_proc_cmd(sample_id,
        //                                   pairs,
        //                                   select_map_type,
        //                                   filter_below_FF,
        //                                   filter_below_RR,
        //                                   filter_below_FR,
        //                                   filter_below_RF,
        //                                   filter_same_frag))
        "touch ${sample_id}.select.pairs.gz"        
}

process filter_stats {
    publishDir "results/filter/stats"
    tag       "${sample_id}"
    container "bskubi/pairtools:1.0.4"
    memory 4.Gb

    input:
        tuple(val(sample_id), path(pairs))
    
    output:
        tuple(val(sample_id), path(pairs), path("${sample_id}.select.stats.txt"))
    
    shell:
        "pairtools stats -o !{sample_id}.select.stats.txt !{pairs}"
}

process cool {
    publishDir "results/cool/full", mode: "copy"
    tag       "${sample_id}"
    memory 4.Gb

    input:
        tuple(val(sample_id), path(pairs), val(assembly), path(chromsizes), val(min_bin_size))
    
    output:
        tuple(val(sample_id), path("${sample_id}.cool"))
    
    shell:
        println("cooler cload pairs --assembly ${assembly} -c1 2 -p1 3 -c2 4 -p2 5 ${chromsizes}:${min_bin_size} ${pairs} ${sample_id}.cool")
        "cooler cload pairs --assembly !{assembly} -c1 2 -p1 3 -c2 4 -p2 5 !{chromsizes}:!{min_bin_size} !{pairs} !{sample_id}.cool"
    
    stub:
        println("cooler cload pairs --assembly ${assembly} -c1 2 -p1 3 -c2 4 -p2 5 ${chromsizes}:${min_bin_size} ${pairs} ${sample_id}.cool")
        "touch ${sample_id}.cool"
}

process downsample {
    publishDir "results/cool/downsample", mode: "copy"
    tag       "${sample_id}"
    memory 12.Gb

    input:
        tuple(val(sample_id), path(pairs), val(count))
    
    output:
        tuple(val(sample_id), path("${sample_id}_${count}.cool"), val("_${count}"))
    
    shell:
        println("cooltools random-sample --exact --count ${count} ${sample_id}.cool ${sample_id}_${count}.cool")
        "cooltools random-sample --exact --count !{count} !{sample_id}.cool !{sample_id}_!{count}.cool"
    
    stub:
        "touch ${sample_id}_${count}.cool"
}

process zoomify {
    publishDir "results/cool", mode: "copy"
    tag       "${sample_id}"
    memory 12.Gb

    input:
        tuple(val(sample_id), path(pairs), val(tag), val(binsizes), path(blacklist))
    
    output:
        tuple(val(sample_id), path("${sample_id}.mcool"), val(tag))
    
    shell:
        "cooler zoomify --nproc 4 --resolutions !{binsizes} --balance --balance-args '--blacklist !{blacklist} --max-iters 1000' --out !{sample_id}.mcool !{sample_id}!{tag}.cool"
    
    stub:
        "touch ${sample_id}_${tag}.mcool"
}

process cooler_merge {
    publishDir "results/cool/merge", mode: "copy"
    tag "${sample_id}"
    memory 12.Gb

    input:
        tuple(val(replicate_id), val(tag), val(pairs))

    output:
        tuple(val(replicate_id), val(tag), path("${replicate_id}${tag}.cool"))
    
    shell:
        filenames = pairs.join(" ")
        "cooler merge !{replicate_id}!{tag}.cool !{filenames}"
    
    stub:
        filenames = pairs.join(" ")
        "touch ${replicate_id}${tag}.cool"
}

process final_zoomify {
    publishDir "results/cool", mode: "copy"
    tag       "${sample_id}"
    memory 12.Gb

    input:
        tuple(val(sample_id), val(tag), path(pairs), val(binsizes), path(blacklist))
    
    output:
        tuple(val(sample_id), path("${sample_id}${tag}.mcool"), val(tag))
    
    shell:
        "cooler zoomify --nproc 4 --resolutions !{binsizes} --balance --balance-args '--blacklist !{blacklist} --max-iters 1000' --out !{sample_id}!{tag}.mcool !{sample_id}!{tag}.cool"
    
    stub:
        println("cooler zoomify --nproc 4 --resolutions ${binsizes} --balance --balance-args '--blacklist ${blacklist} --max-iters 1000' --out ${sample_id}${tag}.mcool ${sample_id}${tag}.cool")
        "touch ${sample_id}${tag}.mcool"
}

workflow {
    channel.fromPath("hisham.csv")
        | splitCsv(header: true, sep: '\t')
        | set{csv}
    
    csv
        | map{format_align_hic_input(it)}
        | set{align_params}

    csv
        | map{tuple(it.sample_id,
                    file(it.chromsizes),
                    null,
                    it.assembly,
                    it.min_mapq,
                    null,
                    it.max_insert_size,
                    null)}
        | set{parse_params}
    
    csv
        | map{tuple(
            it.sample_id,
            ["UU", "UR", "RU"],
            it.filter_below_FF,
            it.filter_below_RR,
            it.filter_below_FR,
            it.filter_below_RF,
            it.filter_same_frag)}
        | set{filter_params}

    csv
        | map{tuple(it.sample_id,
                    it.assembly,
                    file(it.chromsizes),
                    it.min_bin_size)}
        | set{cool_params}
    
    csv
        | map{tuple(it.sample_id,
                    it.downsample)}
        | set{downsample_params}
    
    csv
        | map{tuple(it.sample_id,
                    it.binsizes,
                    file(it.blacklist))}
        | set{mcool_params}
    
    csv
        | map{tuple(it.sample_id,
                    it.replicate)}
        | set{replicate_params}

    align_params
        | align_hic
        | join(parse_params)
        | parse_to_contacts
        | map{tuple(it[0], it[1])}
        | combine(channel.fromPath("resources/hg38_Arima.bed"))
        | map{tuple(it[0], null, it[2], it[1])}
        | restrict
        | join(filter_params)
	| set{filtered}
    
    filtered
        | filter_contacts
        | filter_stats
        | map{tuple(it[0], it[1])}
        | join(cool_params)
        | cool
        | set{cool_full}
    
    cool_full
        | join(downsample_params)
        | downsample
        | set{cool_downsample}
    
    cool_full
        | map{tuple(it[0], it[1], "")}
        | set{cool_full}
    
    cool_full
        | join(mcool_params)
        | set{cool_full_zoomify_input}

    cool_downsample
        | join(mcool_params)
        | concat(cool_full_zoomify_input)
        | set{zoomify_input}

    zoomify_input
        | zoomify
        | set{mcool}
    
    cool_full
        | join(replicate_params)
        | map{tuple(it[-1], it[-2], it[1])}
        | groupTuple(by:[0, 1])
        | set{full_merge}
    
    cool_downsample
        | join(replicate_params)
        | map{tuple(it[-1], it[-2], it[1])}
        | groupTuple(by:[0, 1])
        | set{downsample_merge}
    
    full_merge
        | concat(downsample_merge)
        | cooler_merge
        | set{merged}

    csv
        | map{tuple(it.replicate,
              "",
              it.binsizes,
              file(it.blacklist))}
        | groupTuple(by:[0, 1, 2, 3, 4])
        | set{final_zoomify_params}
    csv
        | map{tuple(it.replicate,
              "_${it.downsample}",
              it.binsizes,
              file(it.blacklist))}
        | groupTuple(by:[0, 1, 2, 3, 4])
        | concat(final_zoomify_params)
        | set{final_zoomify_params}
    
    merged
        | join(final_zoomify_params, by:[0, 1])
        | set{final_zoomify_input}
    
    zoomify_input | view
    final_zoomify_input | view

    final_zoomify_input
        | final_zoomify
}