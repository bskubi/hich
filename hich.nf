workflow keydiff {
    take:
        left
        right
        by
        how

    main:
    /*
    Arguments:
        left -- the left channel treated as a 'table' with LinkedHashMap rows
        right -- the right channel 'table'
        by -- key or list of keys
        how -- 'left' or 'right',
            'left': returns left - right keyset
            'right': returns right - left keyset
        */

    // Start by identifying unique keys from each channel
    // We add true after the key hashmap so that when identifying
    // missing keys with join, we can distinguish matches from missing
    // by having true or null in the second position.
    left_keys = left
        | map{it.subMap(by)}
        | unique
        | map{[it, true]}

    right_keys = right
        | map{it.subMap(by)}
        | unique
        | map{[it, true]}
    
    // Identify which keys in the left are not found in the right
    // Returns [key, true] if a match or [key, null] if no match
    missing =
        left_keys
        | join(right_keys, remainder: true)
        | branch {
            neither: it[1] &&  it[2] //key present in both
            right:   it[1] && !it[2] //key in left only
            left:   !it[1] &&  it[2] //key in right only
        }

    // Get rid of placeholder and keep just the key
    if (how == "right") {
        missing.right | map{it[0]} | set{result}
    } else if (how == "left") {
        missing.left | map{it[0]} | set{result}
    }

    emit:
        result


}

workflow sqljoin {
    take:
    left
    right
    kwargs

    main:
    /* Similar to a SQL inner join

    Tables are represented as channels emitting LinkedHashMap 'rows'

    Arguments:
    left -- left channel
    right -- right channel
    by -- Key or list of keys to join on
    suffix -- suffix to add to right if necessary to avoid non-key clashes

    Notes:
    1. If there is a conflict (i.e. the left table has a non-key 'data'
    and the right has a non-key 'data' and 'data_right' with the suffix
    '_right', then the right keys with higher index have more suffix copies
    appended to them.

    2. sqljoin takes key type into account, so even if key1 == key2,
    they will not be joined (and will not be treated as conflicting and
    thus distinguished with a suffix) unless their class is the same. One
    way this can be an issue is if their types are GString and String.
    */

    by = kwargs.get("by", null)
    suffix = kwargs.get("suffix", "_right")
    how = kwargs.get("how", "left")

    // To accommodate left, right and full joins, we add any missing keys
    if (how == "full" || how == "left") {
        // Add any keys in left that are missing in right to right.
        missing = keydiff(left, right, by, "right")
        right = right | concat(missing)
    }
    if (how == "full" || how == "right") {
        // Add any keys in left that are missing in right to right.
        missing = keydiff(left, right, by, "left")
        left = left | concat(missing)
    }

    // Nextflow cross operator expects the source (left) to have unique
    // keys, so after extracting the keys to match on to a [keys, hashmap]
    // list, we use groupTuple on the left column to get [keys, [hashmaps]]

    left = left | map{[it.subMap(by), it]} | groupTuple
    right = right | map{[it.subMap(by), it]}

    // cross emits one item of format:
    // [[[key, [[leftrow1], [leftrow2], ...]],
    //  [[key, [rightrow]]]
    // for each matching key in the left and right
    joined =
        left
        | cross(right)
        | map{
            result1 = []
            // Drop the key (which is still preserved in the hashmaps)
            // and extract the hashmaps to get:
            // left = [[leftrow1], [leftrow2], ...]
            // rightrow = rightrow
            left = it[0].drop(1)[0]
            rightrow = it[1].drop(1)[0]

            // Iterate through each leftrow and combine with rightrow
            left.each() {
                leftrow ->
                combined = [:]

                // Put all items (including key) from leftrow into combined
                combined.putAll(leftrow)
                
                // Iterate through each key in rightrow,
                // appending suffix if needed to avoid clashes with leftrow
                // and adding it to the combined row.
                rightrow.each() {
                    key, value ->
                    // Don't add keys again as they were added in putAll
                    if (!by.contains(key)) {
                        // Add as many suffixes as necessary to avoid clash
                        while (suffix != "" && combined.containsKey(key)) {
                            key = key + suffix
                        }
                        combined[key] = value
                    }
                }
                result1.add(combined)
            }
            result1
        }
        | collect
        | flatMap{it}
    
    emit:
        joined
}

def getIdx(item, index) {
    if (item instanceof List) {
        validIndex = Math.min(index, item.size() - 1)
        return item[validIndex]
    }
    return item
}

workflow JoinProcessResults {
    take:
        proc
        channels
        input
        output
        join_by
        kwargs
        latest

    main:
        result_jpr = channels[0]
            | map {
                elements = it.subMap(input).values().toList()
                elements = kwargs != false ? elements + it.subMap(kwargs) : elements
                tuple(*elements)
            }
            | proc
            | map{
                result_map = [:];
                [output, it].transpose().each { k, v -> result_map[k] = v};
                latest ? (result_map.latest = result_map[latest]) : null
                result_map
            }


        allchannels = [result_jpr] + channels

        channels.eachWithIndex {ch, i ->
            by = getIdx(join_by, i)
            channels[i] = sqljoin(channels[i], allchannels[i], [by: by, suffix: ""])
            allchannels[i+1] = channels[i]
        }
    
    emit:
        channels[-1]
}

workflow MakeResourceFile {
    take:
        exists
        missing
        no_change
        file_key
        proc
        input
        output
        join_by
        pass_hashmap

    main:
        made = missing
            | map {
                elements = it.subMap(input).values().toList()
                elements = pass_hashmap ? elements + it : elements
                
                tuple(*elements)
            }
            | unique
            | proc
            | map{
                result_jpr_mrf = [:];
                [output, it].transpose().each {k, v -> result_jpr_mrf[k] = v};
                result_jpr_mrf
            }
        
        missing = sqljoin(missing, made, [by: join_by, suffix: ""])

        result_mrf = exists
            | map{
                it[file_key] = file(it[file_key]);
                it
            }
            | concat(no_change)
            | concat(missing)
        
    emit:
        result_mrf
}

process StageReferences {
    input:
    tuple val(assembly), path(url)

    output:
    tuple val(assembly), path(url)

    shell:
    ":"
}

workflow TryDownloadMissingReferences {
    /*
        We will mainly aim to download resources from NCBI:

        https://www.ncbi.nlm.nih.gov/datasets/

        See:
        https://lh3.github.io/2017/11/13/which-human-reference-genome-to-use
    */
    take:
        samples

    main:
        /*  
            1. Convert 'assembly' to url key if no existing reference provided
            2. Reference (URL or local file) cast to file(reference) which
               causes Nextflow to automatically stage URLs for download and
               then treat them as local files going forward.
        */
        synonyms = ["hg38":"hg38",
         "homo_sapiens":"hg38",
         "GRCh38":"hg38",
         "mm10":"mm10",
         "dm6":"dm6",
         "galGal5":"galGal5",
         "bGalGal5":"galGal5",
         "danRer11":"danRer11"]

        /*
            Reference genome URLs for human and for common model organisms
        */
        urls = ["hg38":"https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/@@download/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz",
                "mm10":"https://www.encodeproject.org/files/mm10_no_alt_analysis_set_ENCODE/@@download/mm10_no_alt_analysis_set_ENCODE.fasta.gz",
                "dm6":"https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/215/GCF_000001215.4_Release_6_plus_ISO1_MT/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna.gz",
                "galGal5":"https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/027/408/225/GCA_027408225.1_bGalGal5.pri/GCA_027408225.1_bGalGal5.pri_genomic.fna.gz",
                "danRer11":"https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/035/GCF_000002035.6_GRCz11/GCF_000002035.6_GRCz11_genomic.fna.gz"]
        
        /*
            1. Filter for samples tagged with datatype 'fastq'
            2. Group by 'assembly' and 'reference'
            3. Filter for unspecified or non-existent unique files/URLs
            4. Get URL for missing references
            5. Download
            6. Set file as output path
        */
        
        // Separate to download from those with existing reference
        samples
            | branch{
                download: it.get("reference", "").trim().length() == 0
                          || !file(it.reference).exists()
                exists: true
            }
            | set{branched}
        
        // Convert existing reference to file
        branched.exists
            | map{it.reference = file(it.reference); it}
            | set{exists}
        
        // Map the existing files back to the samples
        
        sqljoin(samples, exists, [by: "sample_id", suffix: ""])
            | set{samples}


        branched.download
            | map{it.subMap("assembly")}
            | map{it.reference = file(urls[synonyms[it.assembly]]); it}
            | unique
            // | StageReferences
            | map{["assembly": it[0], "reference":it[1]]}
            | set{remote_refs}

        download = sqljoin(branched.download, remote_refs, [by: "assembly", suffix: ""])
        
        samples = sqljoin(samples, download, [by: "sample_id", suffix: ""])


    emit:
        samples
}

process MakeChromsizes {
    conda 'bioconda::ucsc-fasize'

    input:
    tuple path(reference), val(assembly), val(chromsizes)

    output:
    tuple val(assembly), path(chromsizes)

    shell:
    "faSize -detailed -tab ${reference} > ${chromsizes}"
}

workflow MakeMissingChromsizes {
    take:
        samples
    
    main:
        def hasChromsizesFilename = {
            it.get("chromsizes", "").toString().trim().length() > 0
            && it.get("chromsizes") != "NULL"}

        def chromsizesExists = {hasChromsizesFilename(it)
                                && file(it.chromsizes).exists()}
        
        def ensureChromsizesFilename = {
            if (!hasChromsizesFilename(it)) {
                it.chromsizes = "${it.assembly}.sizes"
            }
            it
        }

        samples
            | map{ensureChromsizesFilename(it)}
            | branch{
                exists: chromsizesExists(it)
                missing: !chromsizesExists(it) && hasChromsizesFilename(it)
                no_change: true
            }
            | set{branched}
        
        samples = MakeResourceFile(
            branched.exists,
            branched.missing,
            branched.no_change,
            "chromsizes",
            MakeChromsizes,
            ["reference", "assembly", "chromsizes"],
            ["assembly", "chromsizes"],
            ["assembly"],
            false
        )

    emit:
        samples
}

process BwaMem2Index {
    publishDir "resources/index/bwa-mem2/", mode: 'move'
    container "bskubi/bwa-mem2"

    input:
    tuple path(reference), val(prefix)

    output:
    tuple val(prefix), path("${prefix}.0123"), path("${prefix}.amb"),
          path("${prefix}.ann"), path("${prefix}.bwt.2bit.64"),
          path("${prefix}.pac")

    shell:
    "touch ${prefix}.0123 ${prefix}.amb ${prefix}.ann ${prefix}.bwt.2bit.64 ${prefix}.pac"
    //"bwa-mem2 index -p ${prefix} ${reference}"

    stub:
    "touch ${prefix}.0123 ${prefix}.amb ${prefix}.ann ${prefix}.bwt.2bit.64 ${prefix}.pac"
}

workflow MakeMissingIndex {

    take:
        samples
    
    main:

        def fileExists = { dir, file -> new File(dir, file).exists() }
        
        samples
            | filter{it.get("index_dir", "").trim().length() == 0
                     || it.get("index_prefix").trim().length() == 0
                     || !fileExists(it.index_dir, "${it.index_prefix}.0123")}
            | set{no_index}
        
        no_index
            | map{it.subMap("reference", "assembly")}
            | unique
            | map{it.index_prefix = it.index_prefix?.trim() ? it.index_prefix : it.assembly; it}
            | map{tuple(it.reference, it.index_prefix)}
            | BwaMem2Index
            | map{["index_prefix": it[0], "index_dir": "resources/index/bwa-mem2/"]}
            | set {new_index}

        sqljoin(no_index, new_index, [by: "assembly", suffix: ""])
            | set {new_index}

        sqljoin(samples, new_index, [by: "sample_id", suffix: ""])
            | map{it.index_dir = file(it.index_dir); it}
            | set {samples}

    emit:
        samples

}

process MakeDigest {
    input:
    tuple path(reference), val(enzymes), val(fragfile), val(assembly)

    output:
    tuple path(reference), val(enzymes), path(fragfile), val(assembly)

    script:
    "redigest --output ${fragfile} ${reference} ${enzymes}"
}

workflow MakeMissingDigest {
    take:
        samples
    
    main:
        def digestEnzymesDeclared = {it.get("enzymes").trim().length() != 0}

        def hasFragfileName = {it.get("fragfile").trim().length() > 0}
        
        def fragfileExists = {hasFragfileName(it) && file(it.fragfile).exists()}

        def ensureFragfileName = {
            if (digestEnzymesDeclared(it) && !hasFragfileName(it)) {
                it.fragfile = "${it.assembly}_${it.enzymes}.bed"
            }
            it;
        }

        samples
            | map{ensureFragfileName(it)}
            | branch{
                exists: digestEnzymesDeclared(it) && fragfileExists(it)
                missing: digestEnzymesDeclared(it) && !fragfileExists(it)
                no_change: true
            }
            | set{branched}
        
        samples = MakeResourceFile(
            branched.exists,
            branched.missing,
            branched.no_change,
            "fragfile",
            MakeDigest,
            ["reference", "enzymes", "fragfile", "assembly"],
            ["reference", "enzymes", "fragfile", "assembly"],
            ["assembly", "enzymes"],
            false
        )

    emit:
        samples
}

process BwaMem2Align {
    container "bskubi/bwa-mem2"
    //maxRetries 4
    //memory {20.GB + 20.GB * task.attempt}

    // NOTE: Alignment speed is trivially parallelizeable and does not benefit
    // from running alignment in parallel multiple files at once. Each instance
    // of bwa-mem2 uses about 15 gigs of memory. For these two reasons we 
    // tell nextflow to run one alignment process at a time with maxForks 1.
    maxForks 1

    input:
    tuple val(sample_id), path(index_dir), val(index_prefix), path(fastq1), path(fastq2), val(bam)

    output:
    tuple val(sample_id), path(bam)

    shell:
    align = "bwa-mem2 mem -t 10 -SP5M ${index_dir}/${index_prefix} ${fastq1} ${fastq2}"
    tobam = "samtools view -b -o ${bam}"
    "${align} | ${tobam}"
}

workflow Align {
    take:
        samples

    main:
        samples
            | branch {
                fastq: it.datatype == "fastq"
                other: true
            } | set {branched}
        
        branched.fastq
            | map{
                it.data1 = file(it.data1)
                it.data2 = file(it.data2)
                it.bam = "${it.sample_id}.bam"
                it}
            | set {fastq}

        samples = JoinProcessResults(
            BwaMem2Align,
            [fastq, samples],
            ["sample_id", "index_dir", "index_prefix", "data1", "data2", "bam"],
            ["sample_id", "bam"],
            ["sample_id"],
            false,
            "bam")
    emit:
        samples


}

process PairtoolsParse2 {
    container "bskubi/pairtools:1.0.4"

    input:
    tuple val(sample_id), path(bam), path(chromsizes), val(pairs), val(kwargs)

    output:
    tuple val(sample_id), path(pairs)

    shell:
    // Sort by name with sambamba, then parse, then sort by position
    cmd = ["pairtools parse2",
           "--chroms-path ${chromsizes}",
           "--assembly ${kwargs.assembly}",
           "--min-mapq ${kwargs.min_mapq}",
           "--drop-readid",
           "--drop-seq",
           "--drop-sam",
           "${bam}",
           "| pairtools flip --chroms-path ${chromsizes}",
           "| pairtools sort --output ${pairs}"]

    cmd.join(" ")
}

workflow Parse {
    take:
    samples

    main:

        samples
        | map{it.bam = it.datatype in ["sam", "bam"] ? file(it.data1) : file(it.bam); it}
        | filter{"bam" in it && it.bam.exists()}
        | map{it.pairs = "${it.sample_id}.pairs.gz"; it}
        | set {bam}
    
    samples = JoinProcessResults(
        PairtoolsParse2,
        [bam, samples],
        ["sample_id", "bam", "chromsizes", "pairs"],
        ["sample_id", "pairs"],
        ["sample_id"],
        ["assembly", "min_mapq"],
        "pairs")
    
    samples | map{it.id = it.sample_id; it} | set{samples}

    emit:
        samples
}

process Fragtag {

    input:
    tuple val(sample_id), path(pairs), path(fragfile), val(tagged_pairs)

    output:
    tuple val(sample_id), path(tagged_pairs)

    shell:
    // ["pairtools restrict",
    //  "--frags ${fragfile}",
    //  "--output ${tagged_pairs}",
    //  "${pairs}"].join(" ")
    
    cmd = "fragtag ${fragfile} ${tagged_pairs} ${pairs}"
    cmd
}

workflow OptionalFragtag {
    take:
        samples

    main:
        def hasFragfileName = {
            it.get("fragfile").toString().trim().length() > 0
        }
        
        def fragfileExists = {
            hasFragfileName(it) && file(it.fragfile).exists()
        }

        samples
            | filter{fragfileExists(it)}
            | map{it.frag_pairs = "${it.sample_id}_fragtag.pairs.gz"; it}
            | set{fragtag}
        
        samples = JoinProcessResults(
            Fragtag,
            [fragtag, samples],
            ["sample_id", "pairs", "fragfile", "frag_pairs"],
            ["sample_id", "frag_pairs"],
            ["sample_id"],
            false,
            "frag_pairs"
        )
    
    emit:
        samples
}

process Merge {
    //container "bskubi/pairtools:1.0.4"

    input:
    tuple val(id), path(samples)

    output:
    tuple val(id), path("${id}.pairs.gz")

    shell:
    samples = (samples.getClass() == nextflow.processor.TaskPath) ? samples : samples.join(" ")
    "pairtools merge --output ${id}.pairs.gz ${samples}"
}

workflow TechrepsToBioreps {
    take:
        samples

            main:

        def isTechrep = {it.get("is_techrep", "").toString().trim() in ["", "true", true, 1]}
        def hasStructure = {it.get("condition").length() > 0
                            && it.get("biorep").length() > 0}
        def ensureStructure = {
            if (isTechrep(it)) {it.is_techrep = true}
            if (isTechrep(it) && !hasStructure(it)) {
                it.biorep = it.sample_id
                it.condition = it.sample_id
            }
            it
        }

        samples
            | filter{isTechrep(it)}
            | map{ensureStructure(it)}
            | map{tuple(it.subMap("condition", "biorep"), it)}
            | groupTuple
            | map {
                it[0].latest = []
                it[1].each {
                    hashmap ->
                    it[0].latest.add(hashmap.latest)
                }
                it[0].techreps = it[1]
                it[0].is_biorep = true
                it[0].id = "${it[0].condition}_${it[0].biorep}".toString()
                it[0].reference = it[1].reference.unique()[0]
                it[0].chromsizes = it[1].chromsizes.unique()[0]
                it[0]
            }
            | AssignParams
            | set{to_merge}

        to_merge = JoinProcessResults(
            Merge,
            [to_merge],
            ["id", "latest"],
            ["id", "biorep_merge_pairs"],
            ["id"],
            false,
            "biorep_merge_pairs"
        )

        to_merge
            | concat(samples)
            | set {samples}

    emit:
        samples
}

process PairtoolsDedup {
    container "bskubi/pairtools:1.0.4"

    input:
    tuple val(id), path(infile)

    output:
    tuple val(id), path("${id}_dedup.pairs.gz")

    shell:
    //"cp ${infile} ${id}_dedup.pairs.gz "
    "pairtools dedup --output ${id}_dedup.pairs.gz ${infile}"
}

workflow Deduplicate {
    take:
        samples
    
    main:
        
        samples | filter{it.deduplicate} | set {deduplicate}

        JoinProcessResults(
            PairtoolsDedup,
            [deduplicate, samples],
            ["id", "latest"],
            ["id", "dedup_pairs"],
            ["id"],
            false,
            "dedup_pairs"
        ) | set {samples}

    emit:
        samples
}

workflow BiorepsToConditions {
    take:
        samples

    main:

        def isBiorep = {it.containsKey("is_biorep") && it.is_biorep}
        def hasStructure = {it.get("condition").length() > 0
                            && it.get("biorep").length() > 0}
        def ensureStructure = {
            if (isBiorep(it)) {it.is_biorep = true}
            if (isBiorep(it) && !hasStructure(it)) {
                it.biorep = it.sample_id
                it.condition = it.sample_id
            }
            it
        }

        samples
            | filter{isBiorep(it)}
            | map{ensureStructure(it)}
            | map{tuple(it.subMap("condition"), it)}
            | groupTuple
            | map {
                it[0].latest = []
                x = 0
                it[1].each {
                    hashmap ->
                    x += 1
                    it[0].latest += hashmap.latest instanceof ArrayList ? hashmap.latest : [hashmap.latest]
                }
                it[0].reference = it[1].reference.unique()[0]
                it[0].chromsizes = it[1].chromsizes.unique()[0]
                it[0].bioreps = it[1]
                it[0].is_condition = true
                it[0].id = "${it[0].condition}".toString()
                it[0]
            }
            | AssignParams
            | set{to_merge}


        to_merge = JoinProcessResults(
            Merge,
            [to_merge],
            ["id", "latest"],
            ["id", "condition_merge_pairs"],
            ["id"],
            false,
            "condition_merge_pairs"
        )

        to_merge
            | concat(samples)
            | set {samples}

    emit:
        samples
}

process PairtoolsSelect {
    container "bskubi/pairtools:1.0.4"

    input:
    tuple val(id), path(pairs), val(kwargs)

    output:
    tuple val(id), path("${id}_select.pairs.gz")

    shell:

    def quote = { List<String> items ->
    result = items.collect { "'$it'" }.join(", ")
    "[${result}]"

}

    pair_types = "pair_type in ${quote(kwargs.keep_pair_types)}"
    
    cis_trans = kwargs.keep_cis ^ kwargs.keep_trans
                    ? (kwargs.keep_cis
                            ? "(chrom1 == chrom2)"
                            : "(chrom1 != chrom2)")
                    : null

    min_distances = [:]
    min_distances += kwargs.min_dist_fr ? ["+-":kwargs.min_dist_fr] : [:]
    min_distances += kwargs.min_dist_rf ? ["-+":kwargs.min_dist_rf] : [:]
    min_distances += kwargs.min_dist_ff ? ["++":kwargs.min_dist_ff] : [:]
    min_distances += kwargs.min_dist_rr ? ["--":kwargs.min_dist_rf] : [:]
    strand_dist = min_distances.collect {
        strand, dist ->
        s1 = strand[0]
        s2 = strand[1]
        "(strand1 + strand2 == '${s1}${s2}' and abs(pos2 - pos1) >= ${dist})"
    }.join(" or ")
    strand_dist = strand_dist ?: null

    chroms = kwargs.chroms ? "chrom1 in ${quote(kwargs.chroms)} and chrom2 in ${quote(kwargs.chroms)}" : null
    
    frags = kwargs.discard_same_frag ? "rfrag1 != rfrag2" : null
    
    filters = [pair_types, cis_trans, strand_dist, chroms, frags]
    
    filters.removeAll([null])
    filters = filters.collect {"(${it})"}
    filters = filters.join(" and ")

    cmd = """pairtools select --output ${id}_select.pairs.gz "${filters}" ${pairs}"""
    cmd
    
}

workflow Select {
    take:
        samples
    
    main:        
        JoinProcessResults(
            PairtoolsSelect,
            [samples],
            ["id", "latest", "select"],
            ["id", "select_pairs"],
            ["id"],
            false,
            "select_pairs"
        ) | set{samples}


    emit:
        samples
}

process CoolerZoomify {
    conda "cooler"
    maxForks 2

    input:
    tuple val(id), path(infile), path(chromsizes), val(pairs_format), val(matrix)

    output:
    tuple val(id), path("${id}.mcool")

    shell:
    min_bin = matrix.make_matrix_binsizes.min()
    bins = matrix.make_matrix_binsizes.join(",")

    balance = matrix.balance.join(" ")
    balance_args = balance ? "--balance-args '${balance}'" : null
    cmd = ["cooler cload pairs",
           "--chrom1 ${pairs_format.chrom1}",
           "--pos1 ${pairs_format.pos1}",
           "--chrom2 ${pairs_format.chrom2}",
           "--pos2 ${pairs_format.pos2}",
           "${chromsizes}:${min_bin}",
           "${infile} ${id}.cool",
           "&& cooler zoomify",
           "--nproc 10",
           "--resolutions '${bins}'",
           "--balance",
           balance_args,
           "--out ${id}.mcool",
           "${id}.cool"]
    cmd.removeAll([null])
    cmd = cmd.join(" ")
    cmd
}

workflow MakeMcool {
    take:
        samples
    
    main:
        samples
            | filter{it.matrix.make_mcool_file_format}
            | set{mcool}

        JoinProcessResults(
            CoolerZoomify,
            [mcool, samples],
            ["id", "latest", "chromsizes", "pairs_format", "matrix"],
            ["id", "mcool"],
            ["id"],
            false,
            "mcool"
        ) | set{samples}


    emit:
        samples
}

process JuicerToolsPre {
    container "bskubi/juicer_tools:1.22.01"
    maxForks 2

    input:
    tuple val(id), path(infile), path(chromsizes), val(pairs_format), val(matrix)

    output:
    tuple val(id), path("${id}.hic")

    shell:
    outfile = "${id}.hic"
    min_bin = matrix.make_matrix_binsizes.min()
    bins = matrix.make_matrix_binsizes.join(",")

    ["java -Xmx20g -jar /app/juicer_tools.jar pre",
    "${infile} ${outfile} ${chromsizes}"].join(" ")
}

workflow MakeHic {
    take:
        samples
    
    main:
        samples | filter{it.matrix.make_hic_file_format} | set{hic}

        JoinProcessResults(
            JuicerToolsPre,
            [hic, samples],
            ["id", "latest", "chromsizes", "pairs_format", "matrix"],
            ["id", "hic"],
            ["id"],
            false,
            null
        ) | set{samples}

    emit:
        samples
}

workflow AssignParams {
    take:
        samples
    
    main:
        
        samples
            | map {
                sample ->
                params.defaults.each {
                    k, v ->
                    !(k in sample) ? sample += [(k):v] : null
                }
                params.each {
                    k, bundle ->
                    if (bundle.containsKey("ids") && sample.id in bundle.ids) {
                        update = bundle.clone()
                        update.remove("ids")
                        sample += update
                    }
                }
                sample
            }
            | set{samples}
    emit:
        samples
}

workflow {
    channel.fromPath("samples.csv", checkIfExists: true)
        | splitCsv(header: true)
        | map{it.id = it.sample_id; it}
        | AssignParams
        | TryDownloadMissingReferences
        | MakeMissingChromsizes
        | MakeMissingIndex
        | MakeMissingDigest
        | Align
        | Parse
        | OptionalFragtag
        | TechrepsToBioreps
        | Deduplicate
        | BiorepsToConditions
        | Select
        | MakeHic
        | MakeMcool
}

