def shout(message) {
    stars = "*"*message.length()
    // print(stars)
    // print(message)
    // print(stars)
}

def warn(message) {
    print(message)
}

def keydiff (left, right, by, how) {
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
        return missing.right | map{it[0]}
    } else if (how == "left") {
        return missing.left | map{it[0]}
    }
}

def sqljoin (params, left, right) {
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

    by = params.get("by", null)
    suffix = params.get("suffix", "_right")
    how = params.get("how", "left")

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
    
    return joined
}

// def getIdx(item, index) {
//     if (item instanceof List) {
//         validIndex = Math.min(index, item.size() - 1)
//         return item[validIndex]
//     }
//     return item
// }

// def JoinProcessResults(proc, channels, input, output, join_by, kwargs = false, new_latest = null) {
//     if(proc == MergeTechrepsToBioreps) {
//         print("kwargs: ${kwargs}")
//         print("new_latest: ${new_latest}")
//     }
//     shout("JoinProcessResults can't handle processes that don't input or output tuples")
//     result_jpr = channels[0]
//         | map {
//             elements = it.subMap(input).values().toList()
//             elements = kwargs ? elements + it : elements
//             tuple(*elements)
//         }
//         | proc
//         | map{
//             result_map = [:];
//             [output, it].transpose().each { k, v -> result_map[k] = v};
//             result_map
//         }
    
//     if (proc == MergeTechrepsToBioreps) {
//         result_jpr | view
//     }

//     allchannels = [result_jpr] + channels

//     channels.eachWithIndex {ch, i ->
//         by = getIdx(join_by, i)
//         channels[i] = sqljoin(channels[i], allchannels[i], by = by, suffix = "")
//         allchannels[i+1] = channels[i]
//     }
    
//     return channels[-1]
// }

def MakeResourceFile(branched,
                     file_key,
                     proc,
                     input,
                     output,
                     join_by,
                     kwargs = false) {
    shout("MakeResourceFile is untested")

    made = branched.missing
        | map {
            elements = it.subMap(input).values().toList()
            elements = kwargs ? elements + it : elements
            tuple(*elements)
        }
        | unique
        | proc
        | map{
            result_jpr_mrf = [:];
            [output, it].transpose().each {k, v -> result_jpr_mrf[k] = v};
            result_jpr_mrf
        }
    
    missing = sqljoin(branched.missing, made, by: join_by, suffix: "")

    result_mfr = branched.exists
        | map{
            it[file_key] = file(it[file_key]);
            it
        }
        | concat(branched.no_change)
        | concat(missing)
    
    return result_mfr
}

// def GroupHashmaps(hashmaps) {
//     grouped = [:]
//     hashmaps.each {
//         hashmap ->
//         hashmap.keySet().each {
//             k ->
//             grouped[k] = []
//         }
//     }
//     hashmaps.each {
//         hashmap ->
//         hashmap.each {
//             k, v ->
//             grouped[k].add(v)
//         }
//     }
//     grouped.each {
//         k, v ->
//         grouped[k] = grouped[k].unique()
//         if (grouped[k].size() == 1) {
//             grouped[k] = grouped[k][0]
//         }
//     }

//     return grouped
// }

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
        shout("WARNING: StageReferences is commented out in EnsureReferences")
        shout("WARNING: Have not tested ability to recombine branched samples in EnsureReferences")
        
        samples
            | branch{
                download: it.get("reference", "").trim().length() == 0
                          || !file(it.reference).exists()
                exists: true
            }
            | set{branched}
        
        branched.exists
            | map{it.reference = file(it.reference); it}
            | set{exists}
        
        samples = sqljoin(exists, samples, by: "sample_id", suffix: "")

        
        branched.download
            | map{it.subMap("assembly")}
            | map{it.reference = file(urls[synonyms[it.assembly]]); it}
            | unique
            // | StageReferences
            | map{["assembly": it[0], "reference":it[1]]}
            | set{remote_refs}

        download = sqljoin(branched.download, remote_refs, by: "assembly", suffix: "")
        
        samples = sqljoin(samples, download, by: "sample_id", suffix: "")
        
        

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
                missing: hasChromsizesFilename(it)
                no_change: true
            }
            | set{branched}
        
        samples = MakeResourceFile(
            branched,
            "chromsizes",
            MakeChromsizes,
            ["reference", "assembly", "chromsizes"],
            ["assembly", "chromsizes"],
            ["assembly"]
        )

    emit:
        samples
}

// process MakeDigest {
//     input:
//     tuple path(reference), val(enzymes), val(fragfile), val(assembly)

//     output:
//     tuple path(reference), val(enzymes), path(fragfile), val(assembly)

//     script:
//     "redigest --output ${fragfile} ${reference} ${enzymes}"
// }

// workflow MakeMissingDigest {
//     take:
//         samples
    
//     main:
//         def digestEnzymesDeclared = {it.get("enzymes").trim().length() != 0}

//         def hasFragfileName = {it.get("fragfile").trim().length() > 0}
        
//         def fragfileExists = {hasFragfileName(it) && file(it.fragfile).exists()}

//         def ensureFragfileName = {
//             if (digestEnzymesDeclared(it) && !hasFragfileName(it)) {
//                 it.fragfile = "${it.assembly}_${it.enzymes}.bed"
//             }
//             it;
//         }

//         samples
//             | map{ensureFragfileName(it)}
//             | branch{
//                 exists: digestEnzymesDeclared(it) && fragfileExists(it)
//                 missing: digestEnzymesDeclared(it) && !fragfileExists(it)
//                 no_change: true
//             }
//             | set{branched}
        
//         samples = MakeResourceFile(
//             branched,
//             "fragfile",
//             MakeDigest,
//             ["reference", "enzymes", "fragfile", "assembly"],
//             ["reference", "enzymes", "fragfile", "assembly"],
//             ["assembly", "enzymes"]
//         )

//     emit:
//         samples
// }

// process BwaMem2Index {
//     publishDir "resources/index/bwa-mem2/", mode: 'move'
//     container "bskubi/bwa-mem2"

//     input:
//     tuple path(reference), val(prefix)

//     output:
//     tuple val(prefix), path("${prefix}.0123"), path("${prefix}.amb"),
//           path("${prefix}.ann"), path("${prefix}.bwt.2bit.64"),
//           path("${prefix}.pac")

//     shell:
//     "touch ${prefix}.0123 ${prefix}.amb ${prefix}.ann ${prefix}.bwt.2bit.64 ${prefix}.pac"
//     //"bwa-mem2 index -p ${prefix} ${reference}"

//     stub:
//     "touch ${prefix}.0123 ${prefix}.amb ${prefix}.ann ${prefix}.bwt.2bit.64 ${prefix}.pac"
// }

// workflow MakeMissingIndex {

//     take:
//         samples
    
//     main:
        

//         shout("BWA-MEM2 indexing is commented out in 'shell' section")
//         shout("Indexing only works for BWA-MEM2")

//         def fileExists = { dir, file -> new File(dir, file).exists() }
        
//         samples
//             | filter{it.get("index_dir", "").trim().length() == 0
//                      || it.get("index_prefix").trim().length() == 0
//                      || !fileExists(it.index_dir, "${it.index_prefix}.0123")}
//             | set{no_index}
        
//         no_index
//             | map{it.subMap("reference", "assembly")}
//             | unique
//             | map{it.index_prefix = it.index_prefix?.trim() ? it.index_prefix : it.assembly; it}
//             | map{tuple(it.reference, it.index_prefix)}
//             | BwaMem2Index
//             | map{["index_prefix": it[0], "index_dir": "resources/index/bwa-mem2/"]}
//             | set {new_index}

//         sqljoin(no_index, new_index, by = "assembly", suffix = "")
//             | set {new_index}

//         sqljoin(samples, new_index, by = "sample_id", suffix = "")
//             | map{it.index_dir = file(it.index_dir); it}
//             | set {samples}

//     emit:
//         samples

// }

// process BwaMem2Align {
//     container "bskubi/bwa-mem2"
//     //maxRetries 4
//     //memory {20.GB + 20.GB * task.attempt}

//     // NOTE: Alignment speed is trivially parallelizeable and does not benefit
//     // from running alignment in parallel multiple files at once. Each instance
//     // of bwa-mem2 uses about 15 gigs of memory. For these two reasons we 
//     // tell nextflow to run one alignment process at a time with maxForks 1.
//     maxForks 1

//     input:
//     tuple val(sample_id), path(index_dir), val(index_prefix), path(fastq1), path(fastq2), val(bam)

//     output:
//     tuple val(sample_id), path(bam)

//     shell:
//     align = "bwa-mem2 mem -t 10 -SP5M ${index_dir}/${index_prefix} ${fastq1} ${fastq2}"
//     tobam = "samtools view -b -o ${bam}"
//     "${align} | ${tobam}"
// }

// workflow Align {
//     take:
//         samples

//     main:
//         shout("Only aligns with BWA-MEM2")
//         shout("Add bam squeeze to alignment step")
//         samples
//             | branch {
//                 fastq: it.datatype == "fastq"
//                 other: true
//             } | set {branched}
        
//         branched.fastq
//             | map{
//                 it.data1 = file(it.data1)
//                 it.data2 = file(it.data2)
//                 it.bam = "${it.sample_id}.bam"
//                 it}
//             | set {fastq}

//         samples = JoinProcessResults(BwaMem2Align,
//             [fastq, samples],
//             input = ["sample_id", "index_dir", "index_prefix", "data1", "data2", "bam"],
//             output = ["sample_id", "bam"],
//             join_by = ["sample_id"])
//     emit:
//         samples


// }

// process PairtoolsParse2 {
//     container "bskubi/pairtools:1.0.4"

//     input:
//     tuple val(sample_id), path(bam), path(chromsizes), val(pairs), val(kwargs)

//     output:
//     tuple val(sample_id), path(pairs)

//     shell:
//     shout("Numerous flags currently not settable for pairtools parse2")
//     shout("Make sambamba + pairtools dockerfile")
//     // Sort by name with sambamba, then parse, then sort by position
//     cmd = ["pairtools parse2",
//            "--chroms-path ${chromsizes}",
//            "--assembly ${kwargs.assembly}",
//            "--min-mapq 30",
//            "--drop-readid",
//            "--drop-seq",
//            "--drop-sam",
//            "${bam}",
//            "| pairtools sort",
//            "--output ${pairs}"]
//     cmd.join(" ")
// }

// workflow Parse {
//     take:
//     samples

//     main:
//     shout("PairtoolsParse2 doesn't sort by name")

//     samples
//         | map{it.bam = it.datatype in ["sam", "bam"] ? file(it.data1) : file(it.bam); it}
//         | filter{"bam" in it && it.bam.exists()}
//         | map{it.pairs = "${it.sample_id}.pairs.gz"; it}
//         | set {bam}
    
//     samples = JoinProcessResults(
//         PairtoolsParse2,
//         [bam, samples],
//         input = ["sample_id", "bam", "chromsizes", "pairs"],
//         output = ["sample_id", "pairs"],
//         join_by = ["sample_id"],
//         kwargs = true)
    
//     emit:
//     samples
// }

// process Fragtag {
//     container "bskubi/pairtools:1.0.4"
//     maxForks 1

//     input:
//     tuple val(sample_id), path(pairs), path(fragfile), val(tagged_pairs)

//     output:
//     tuple val(sample_id), path(tagged_pairs)

//     shell:
//     ["pairtools restrict",
//      "--frags ${fragfile}",
//      "--output ${tagged_pairs}",
//      "${pairs}"].join(" ")
// }

// workflow OptionalFragtag {
//     take:
//         samples

//     main:
//         def hasFragfileName = {
//             it.get("fragfile").toString().trim().length() > 0
//         }
        
//         def fragfileExists = {
//             hasFragfileName(it) && file(it.fragfile).exists()
//         }

//         samples
//             | filter{fragfileExists(it)}
//             | map{it.frag_pairs = "${it.sample_id}_fragtag.pairs.gz"; it}
//             | set{fragtag}
        
//         samples = JoinProcessResults(
//         Fragtag,
//         [fragtag, samples],
//         input = ["sample_id", "pairs", "fragfile", "frag_pairs"],
//         output = ["sample_id", "frag_pairs"],
//         join_by = ["sample_id"],
//         kwargs = false)
    
//     emit:
//         samples
// }

// process MergeTechrepsToBioreps {
//     container "bskubi/pairtools:1.0.4"

//     input:
//     tuple val(condition), val(biorep), path(samples)

//     output:
//     tuple val(condition), val(biorep), path("${condition}_${biorep}.pairs.gz")

//     shell:
//     "pairtools merge --output ${condition}_${biorep}.pairs.gz ${samples.join(" ")}"
// }

// workflow TechrepsToBioreps {
//     take:
//         samples

//     main:

//         def isTechrep = {it.get("is_techrep", "").toString().trim() in ["", "true", true, 1]}
//         def hasStructure = {it.get("condition").length() > 0
//                             && it.get("biorep").length() > 0}
//         def ensureStructure = {
//             if (isTechrep(it)) {it.is_techrep = true}
//             if (isTechrep(it) && !hasStructure(it)) {
//                 msg = [
//                     "${it.sample_id} listed as techrep, but no biorep and condition given.",
//                     "Treating as a unique condition with biorep=${sample_id} and condition=${sample_id}."
//                 ].join(" ")
//                 warn(msg)
//                 it.biorep = it.sample_id
//                 it.condition = it.sample_id
//             }
//             it
//         }

//         samples
//             | filter{isTechrep(it)}
//             | map{ensureStructure(it)}
//             | map{tuple(it.subMap("condition", "biorep"), it)}
//             | groupTuple
//             | map {
//                 it[0].latest = []
//                 it[1].each {
//                     hashmap ->
//                     keys = ["frag_pairs", "pairs"]
//                     for (k in keys) {
//                         if (k in hashmap.keySet()) {
//                             it[0].latest.add(hashmap[k])
//                             break
//                         }
//                     }
//                 }
//                 it[0].techreps = it[1]
//                 it[0].is_biorep = true
//                 it[0]
//             }
//             | set{to_merge}


//         to_merge = JoinProcessResults(
//         MergeTechrepsToBioreps,
//         [to_merge],
//         input = ["condition", "biorep", "latest"],
//         output = ["condition", "biorep", "biorep_merge_pairs"],
//         join_by = [["condition", "biorep"]],
//         new_latest = "biorep_merge_pairs")

//         to_merge
//             | map{it.latest = it.biorep_merge_pairs; it}
//             | set{to_merge}

//         samples | concat(to_merge) | set{samples}

//     emit:
//         samples
// }

// process Deduplicate {
//     container "bskubi/pairtools:1.0.4"

//     input:
//     tuple path(infile), val(outfile), val(data)

//     output:
//     tuple path(outfile), val(data)

//     shell:
//     "pairtools dedup --output ${outfile} ${infile}"
// }

// workflow DeduplicateBioreps {
//     take:
//         samples
    
//     main:
//         samples
//             | map{it.latest}
//             | view

//     emit:
//         samples
// }

workflow {

    channel.fromPath("samples.csv", checkIfExists: true)
        | splitCsv(header: true)
        | TryDownloadMissingReferences
        | MakeMissingChromsizes
        // | MakeMissingIndex
        // | MakeMissingDigest
        // | Align
        // | Parse
        // | OptionalFragtag
        // | TechrepsToBioreps
        // | DeduplicateBioreps

}

