
workflow keydiff {
    take:
        left
        right
        by
        how

    main:
    /*
        Motivation and implementation logic:

        When Nextflow channel items are LinkedHashMaps ("hashmaps"), they are
        similar to un-normalized tables. Items are rows, keys represent the
        column name for that row, and values represent cell entries.

        We can therefore do an SQL-like join on them. A left join, for example,
        keeps all 'left' rows, adding 'right' rows where there is a matching]
        key and null or default values where there is no matching 'right' key.

        This workflow is used to determine which key column values are missing
        in the left table, right table, or in neither table.

        'by': the set of hashmap keys (i.e. table columns) to join on
        'keyset': the values in the by keys of a particular channel item (i.e. table row)
        
    */

    /* Extract unique keysets from left and right tables.

    This also adds a 'true' value which will be used in the next step to detect
    unmatched keysets.
    */
    left_keys = left
        | map{it.subMap(by)}
        | unique
        | map{[it, true]}

    right_keys = right
        | map{it.subMap(by)}
        | unique
        | map{[it, true]}
    
    /* Determine left ∩ right, left - right, and right - left keysets
    
    Nextflow join operator matches on first element of each item (the keyset)
    The remainder argument keeps unmatched items, regardless of which
    channel they were in. It adds a null value as a placeholder which is groovy
    falsey. If the keyset is missing in neither (found in both), the output
    will be [keyset, true, true] from the true values added in the previous
    step. If the keyset is missing in the left, [keyset, null, true] will be
    emitted, or [keyset, true, null] if it is missing in the right.
    */
    missing =
        left_keys
        | join(right_keys, remainder: true)
        | branch {
            neither: it[1] &&  it[2] //key present in both
            right:   it[1] && !it[2] //key in left only
            left:   !it[1] &&  it[2] //key in right only
        }

    /*
        Extract the keyset for rows found only in the left or only in the right
        depending on the 'how' argument, so that default values can be added
        where necessary during the join.
    */
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
    condition

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

    /*
        Motivation and implementation logic:

        When Nextflow channel items are LinkedHashMaps ("hashmaps"), they are
        similar to un-normalized tables. Items are rows, keys represent the
        column name for that row, and values represent cell entries.

        We can therefore do an SQL-like join on them. A left join, for example,
        keeps all 'left' rows, adding 'right' rows where there is a matching]
        key and null or default values where there is no matching 'right' key.
        If two or more 'right' rows match one 'left' row, then multiple
        rows are produced containing one copy of the 'left' items and one
        copy of one of the matching 'right' rows' items.

        Nextflow's 'join' is similar but importantly different from an SQL
        inner join. Other join types are not available. None operate on named
        keys from hashmaps. We implement a closer approximation to the SQL
        full, left, right and inner joins on Nextflow channels containing
        a single LinkedHashMap item here.

        'keyset': the values in the by keys of a particular channel item (i.e. table row)

        The implementation idea is as follows:

        1. Depending on join type, determine missing keysets that would
        inappropriately lead to rows being dropped and add in empty placeholders.
        For example, for a left join all left rows should be kept. So any keysets
        in the left table not present in the right table must be added as
        placeholders to the right table.

        2. Group the left channel by keyset, giving a [keyset: [list of items]]
        channel. This is necessary because the Nextflow join operator returns
        only the first item from the left channel for each specific keyset.

        3. Use Nextflow cross operator to get pairwise combinations of left
        channel items and right channel items with matching keysets. For example,
        if there are two left items and two right items with a particular
        keyset we obtain:
            [[keyset, [left_1, left_2]], [keyset, right_1]]
            [[keyset, [left_1, left_2]], [keyset, right_2]]

        4. Reformat to produce hashmap-formatted elements:
            [left_1 + right_1]
            [left_1 + right_1]
            [left_1 + right_1]
            [left_1 + right_1]

            Note: the matching keyset will be represented once in the
            resulting hashmaps

            Where there are non-keyset matching columsn in the left and right
            tables, we add a suffix to avoid a clash or (if a blank suffix
            is provided) we overwrite.
    */

    by = condition.get("by", null)
    suffix = condition.get("suffix", "_right")
    how = condition.get("how", "left")

    /* 1. Depending on join type, determine missing keysets that would
       inappropriately lead to rows being dropped and add in empty placeholders.
       For example, for a left join all left rows should be kept. So any keysets
       in the left table not present in the right table must be added as
       placeholders to the right table.
    */
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

    /* 2. Group the left channel by keyset, giving a [keyset: [list of items]]
       channel. This is necessary because the Nextflow join operator returns
       only the first item from the left channel for each specific keyset.
    */

    left = left | map{[it.subMap(by), it]} | groupTuple
    right = right | map{[it.subMap(by), it]}

    /* 3. Use Nextflow cross operator to get pairwise combinations of left
       channel items and right channel items with matching keysets. For example,
       if there are two left items and two right items with a particular
       keyset we obtain:
           [[keyset, [left_1, left_2]], [keyset, right_1]]
           [[keyset, [left_1, left_2]], [keyset, right_2]]
    */

    left
    | cross(right)
    | map{
        /* Reformat item format
           
           from [[keyset, [left_1, left_2]], [keyset, right]]

           to   [left_1 + right]
                [left_2 + right]

            These are hashmaps with keyset represented only once
        */
        result = []

        /* Drop the key (which is still preserved in the hashmaps)
           and extract the hashmaps from the [keyset, hashmap] elements to get:
           left = [leftrow1, leftrow2, ...]
           rightrow = rightrow

           it[0] = [keyset, [left_1, left_2, ...]]
           it[1] = [keyset, right]
           
           drop(1) means to drop 1 item from the beginning of the list, so it
           will drop the keyset.
        */

        /* This is an old implementation. I'm pretty sure it is unnecessary
        to drop and then extract the first element, but I'm leaving it in
        case it introduces bugs.
        */
        //left = it[0].drop(1)[0]       
        //rightrow = it[1].drop(1)[0]
        
        /* The new version just extracts the desired elements directly*/
        left = it[0][1]
        rightrow = it[1][1]

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

                // Add non-keyset keys from right to left.
                if (!by.contains(key)) {

                    // Add as many suffixes as necessary to avoid clash
                    while (suffix != "" && combined.containsKey(key)) {
                        key = key + suffix
                    }

                    /* Add the right key value if either criterion holds:
                        - There is no name clash, or
                        - The left table is not declared dominant
                    */
                    if (condition.get("dominant") != "left" || !combined.containsKey(key)) {
                        combined[key] = value
                    }
                }
            }
            // Add the new hashmap to the list of hashmaps
            result.add(combined)
        }

        // List of per-keyset joins
        result
    }
    | collect       // Currently, channel items are lists of per-keyset joins.
    | flatMap{it}   // Collects to a single list of lists, then flattens
    | set{joined}   // Yielding the desired result of single-hashmap join items.
    
    emit:
        joined
}

def getIdx(item, index) {
    /* If item is a list, get value at index, where too-large indexes are set
       to index of final item. If item is not a list, return it.
    */
    if (item instanceof List) {
        validIndex = Math.min(index, item.size() - 1)
        return item[validIndex]
    }
    return item
}

workflow JoinProcessResults {
    /*
        proc -- the process to call
        channels -- a list of channels.
            The first channel contains single-hashmap items that have as a
            subset of their keys the required process inputs.

            Outputs from the process will be rejoined to channel[0], then
            channel[0] will be joined to channel[1], channel[1] to channel[1], ...
            and the result of the final join will be emitted.
        input -- a list of keys to extract from channel[0] hashmap items and
            use as inputs to the process
        output -- a list of keys in the order expected from the process
        join_by -- either a non-List universal join key or a List of join keys
               Iterate through the channels pairwise ([ch0, ch1], [ch1, ch2]...)
               Then join the first channel to the second. The caller has several
               ways to pass join by parameters:
                    If a List is passed as the join by parameter, then the ith
                    pair of channels are joined on the ith element or on the final
                    join by parameter, whichever is lower.

                    If a non-List is passed, then every pair of channels are
                    joined by it.

                    If it's desired to use a single List L as the join by parameter
                    for every pair of channels, simply pass [L] as the join by
                    parameter (i.e. a single-element list containing only L) 
        condition -- if non-falsey, used to subMap channel[0] and passed as a val
                  parameter to the process
        latest -- if non-falsey, sets the 'latest' parameter of the channel items
                  used as inputs to the process to the specified output
                  from the process.

            Motivation:

            Using hashmaps to keep track of process inputs and outputs has
            numerous advantages, but adds complications.
            
            First, although Nextflow can accept hashmaps as val() params,
            it won't stage files stored in them to the work directory, making
            those files inaccessible to the process.

            Second, any changes to process inputs trigger a rerun of the process.
            We want to avoid this in cases where changes are made to hashmap
            keys that are not used by the process.

            This requires extracting the specific desired elements and passing
            them to the process, but this in turn means that the process can't
            just emit the full hashmap. We therefore have to have a way of
            extracting the desired elements, passing them to the process,
            and rejoining the process outputs to the original hashmap. This
            is what the JoinProcessResults workflow accomplishes.
    */
    take:
        proc
        channels
        input
        output
        join_by
        condition
        latest

    main:
        result = channels[0]
            | map {
                /* Extract needed process inputs in order. Sometimes it is
                convenient to pass an additional set of arguments as a final
                value parameter, which is what the second line does. Then
                format them as a tuple and pass to the process.
                */
                elements = it.subMap(input).values().toList()
                elements = condition ? elements + it.subMap(condition) : elements
                tuple(*elements)
            }
            | proc  // Call the process
            | map{
                /* Collect the outputs from the process and put them into
                a hashmap, using the order specified in 'output'. The transpose
                operator is similar to Python's zip operator, generating
                [[output[0], proc_outputs[0]], [output[1], proc_outputs[1]], ...]
                which are then combined into a hashmap [output[0]:proc_outputs[0], ...]
                during the each{} loop.

                If a non-falsey latest parameter is specified, we set one of the
                process outputs as the value of the latest parameter. This helps
                keep track of the most recent output to facilitate things like
                certain samples skipping processing steps.
                */
                proc_outputs ->
                result_map = [:];
                [output, proc_outputs].transpose().each { k, v -> result_map[k] = v};
                latest ? (result_map.latest = result_map[latest]) : null
                result_map
            }

        /*
            1. Create an extended list of channels, starting wtih the output
                hashmaps from the process.
            
            2. Iterate through the channels pairwise ([ch0, ch1], [ch1, ch2]...)
               Then join the first channel to the second. The caller has several
               ways to pass join by parameters:
                    If a List is passed as the join by parameter, then the ith
                    pair of channels are joined on the ith element or on the final
                    join by parameter, whichever is lower.

                    If a non-List is passed, then every pair of channels are
                    joined by it.

                    If it's desired to use a single List L as the join by parameter
                    for every pair of channels, simply pass [L] as the join by
                    parameter (i.e. a single-element list containing only L) 
        */
        allchannels = [result] + channels

        channels.eachWithIndex {ch, i ->
            by = getIdx(join_by, i)
            channels[i] = sqljoin(channels[i], allchannels[i], [by: by, suffix: ""])
            allchannels[i+1] = channels[i]
        }
    
    emit:
        // The final channel has the complete join specified by the user.
        channels[-1]
}

workflow MakeResourceFile {
    /*
        Motivation:
        
        There are several types of common resource files we may wish to download
        or produce from another resource file. Examples include reference genomes,
        chromsizes files, and restriction digest .bed files. Often multiple
        samples will use common reference files. This workflow identifies the
        set of common references needed by all the samples, uses the given process
        to obtain them, and then adds the file to all the samples that need it
        to avoid redundant downloading and processing.

        If the common resource files are missing, they are created. If they
        exist no special processing is done but it is inforced that they are
        of type 'file' rather than type 'string'. If no_change is required they
        are left as-is.

        exists -- channel of resource files specs that exist and where the values
        under 'file_key' should be converted to type 'file'
        missing -- channel of resource file specs that do not exist and need to be made
        no_change -- channel of resource files that should not be altered
        file_key -- list of keys where the subMap of hashmaps in 'exists' with these
        keys will be converted to type 'file'
        proc -- the Nextflow process to run to make missing files
        input -- the list of inputs expected by the process
        output -- the names of the outputs produced by the process
        join_by -- keys to join the 'made' outputs to the 'missing' items
    */
    take:
        exists
        missing
        no_change
        file_key
        proc
        input
        output
        join_by

    main:
        missing
        | map {
            // Extract input value sfrom hashmap in order required by process
            elements = it.subMap(input).values().toList()
            tuple(*elements)
        }
        | unique    // Get unique inputs to avoid redundant processing
        | proc      // Call process
        | map{
            // Use given output names and process outputs as key:value pairs
            // in a hashmap
            proc_outputs ->
            result = [:]
            [output, proc_outputs].transpose().each {k, v -> result[k] = v}
            result
        }
        | set{made}
        
        // Join the process outputs
        missing = sqljoin(missing, made, [by: join_by, suffix: ""])

        // Make files specified by 'file_key' into files in the 'exists' channel
        // then concatenate with the made and unchanged resource files.
        result = exists
            | map{
                
                it[file_key] = file(it[file_key]);
                it
            }
            | concat(no_change)
            | concat(missing)
        
    emit:
        result
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

process PairtoolsStats {
    publishDir params.general.publish.pair_stats ? params.general.publish.pair_stats : "results",
               saveAs: {params.general.publish.pair_stats ? it : null}
    container "bskubi/pairtools:1.0.4"

    input:
    tuple val(id), val(pairs_id), path(pairs)

    output:
    tuple val(id), val(pairs_id), path("${id}.${pairs_id}.stats.txt")

    shell:
    "pairtools stats --output ${id}.${pairs_id}.stats.txt ${pairs}"
}

process MultiQC {
    publishDir params.general.publish.qc ? params.general.publish.qc : "results",
               saveAs: {params.general.publish.qc ? it : null}
    conda "multiqc.yaml"

    input:
    tuple val(report_name), path(stats)

    output:
    path("${report_name}.multiqc_report.html")

    shell:
    "multiqc --force --filename ${report_name}.multiqc_report.html --module pairtools ."
}

workflow QCReads {
    take:
    samples
    report_name
    
    main:
    samples
    | map{
        sample ->
        pairs = ["pairs",
                  "frag_pairs",
                  "dedup_pairs",
                  "select_pairs",
                  "biorep_merge_pairs",
                  "condition_merge_pairs"]
        
        emit = pairs.collect{
            pairfile ->
            sample.containsKey(pairfile) ? [sample.id, pairfile, sample[pairfile]] : null
        }.findAll{it != null}
        emit
    }
    | collect
    | flatMap
    | PairtoolsStats
    | map{it[2]}
    | collect
    | map{[report_name, it]}
    | MultiQC
    

    emit:
    samples

}

process StageReferences {
    /*  When a URL is passed to a Nextflow function, the resource will be
        automatically downloaded and staged by Nextflow. This is a dummy
        function used to download a unique reference and intentionally has
        no content.
    */
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
            | StageReferences
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
            ["assembly"]
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
            ["assembly", "enzymes"]
        )

    emit:
        samples
}

process BwaMem2Align {
    publishDir params.general.publish.bam ? params.general.publish.bam : "results",
               saveAs: {params.general.publish.bam ? it : null}

    container "bskubi/bwa-mem2"
    //maxRetries 4
    //memory {20.GB + 20.GB * task.attempt}

    // NOTE: Alignment speed is trivially parallelizeable and does not benefit
    // from running alignment in parallel multiple files at once. Each instance
    // of bwa-mem2 uses about 15 gigs of memory. For these two reasons we 
    // tell nextflow to run one alignment process at a time with maxForks 1.
    maxForks 1

    input:
    tuple val(sample_id), path(index_dir), val(index_prefix), path(fastq1), path(fastq2)

    output:
    tuple val(sample_id), path("${sample_id}.bam")

    shell:
    align = "bwa-mem2 mem -t 10 -SP5M ${index_dir}/${index_prefix} ${fastq1} ${fastq2}"
    tobam = "samtools view -b -o ${sample_id}.bam"
    "${align} | ${tobam}"
}

workflow Align {
    take:
        samples

    main:
       
        samples
            | filter{it.datatype == "fastq"
                     && it.get("fastq1")
                     && it.get("fastq2")
            }
            | map{
                it.fastq1 = file(it.fastq1)
                it.fastq2 = file(it.fastq2)
                it
            }
            | set {fastq}

        samples = JoinProcessResults(
            BwaMem2Align,
            [fastq, samples],
            ["sample_id", "index_dir", "index_prefix", "fastq1", "fastq2"],
            ["sample_id", "sambam"],
            ["sample_id"],
            false,
            "sambam")

        if (params.general.get("last_step") == "align") {
            channel.empty() | set{samples}
        }
    emit:
        samples


}

process PairtoolsParse2 {
    publishDir params.general.publish.parse ? params.general.publish.parse : "results",
               saveAs: {params.general.publish.parse ? it : null}
    container "bskubi/pairtools:1.0.4"

    input:
    tuple val(sample_id), path(sambam), path(chromsizes), val(assembly), val(parse_params)

    output:
    tuple val(sample_id), path("${sample_id}.pairs.gz")

    shell:
    cmd = ["pairtools parse2",
           "--assembly ${assembly}",
           "--chroms-path ${chromsizes}"] +
           parse_params +
          ["${sambam} | pairtools sort --output ${sample_id}.pairs.gz"]
    cmd.removeAll([null])

    cmd.join(" ")
}

workflow Parse {
    take:
    samples

    main:

        samples
        | filter{it.get("sambam") && file(it.sambam).exists()}
        | map{it.sambam = file(it.sambam); it}
        | set {sambam}
    
    samples = JoinProcessResults(
        PairtoolsParse2,
        [sambam, samples],
        ["sample_id", "sambam", "chromsizes", "assembly", "parse_params"],
        ["sample_id", "pairs"],
        ["sample_id"],
        null,
        "pairs")
    
    samples | map{it.id = it.sample_id; it} | set{samples}

    if ("Parse" in params.general.get("qc_after")) {
        QCReads(samples, "Parse")
    }

    if (params.general.get("last_step") == "parse") {
        channel.empty() | set{samples}
    }

    emit:
        samples
}

process PairtoolsFlipSort {
    publishDir params.general.publish.flip_sort ? params.general.publish.flip_sort : "results",
               saveAs: {params.general.publish.flip_sort ? it : null}
    container "bskubi/pairtools:1.0.4"

    input:
    tuple val(sample_id), path(pairs), path(chromsizes)

    output:
    tuple val(sample_id), path("${sample_id}.pairs.gz")

    shell:
    cmd = ["pairtools flip --chroms-path ${chromsizes} ${pairs}",
           "| pairtools sort --output ${sample_id}.pairs.gz"]
    cmd.removeAll([null])
    cmd.join(" ")
}

workflow IngestPairs {
    take:
        samples

    main:
        samples
            | filter{it.datatype == "pairs" && it.get("pairs") && file(it.get("pairs")).exists()}
            | map{
                it.pairs = file(it.pairs);
                it.id = it.sample_id;
                it
            }
            | set{ingest}
        
        samples = JoinProcessResults(
            PairtoolsFlipSort,
            [ingest, samples],
            ["sample_id", "pairs", "chromsizes"],
            ["sample_id", "pairs"],
            ["sample_id"],
            null,
            "pairs")
        
        if (params.general.get("last_step") == "IngestPairs") {
            channel.empty() | set{samples}
        }

        if ("IngestPairs" in params.general.get("qc_after")) {
            QCReads(samples, "IngestPairs")
        }

    emit:
        samples
}

process Fragtag {
    publishDir params.general.publish.fragtag ? params.general.publish.fragtag : "results",
               saveAs: {params.general.publish.fragtag ? it : null}


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

        if ("OptionalFragtag" in params.general.get("qc_after")) {
            QCReads(samples, "OptionalFragtag")
        }

    emit:
        samples
}

process Merge {
    publishDir params.general.publish.fragtag ? params.general.publish.fragtag : "results",
               saveAs: {params.general.publish.fragtag ? it : null}

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
                biorep, techreps ->
                biorep += techreps[0].findAll{entry -> techreps.every {map -> map[entry.key] == entry.value}}
                biorep.latest = techreps.collect{it.latest}
                biorep.is_biorep = true
                biorep.id = "${biorep.condition}_${biorep.biorep}"
                biorep
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

        if ("TechrepsToBioreps" in params.general.get("qc_after")) {
            QCReads(samples, "TechrepsToBioreps")
        }

    emit:
        samples
}

process PairtoolsDedup {
    publishDir params.general.publish.dedup ? params.general.publish.dedup : "results",
               saveAs: {params.general.publish.dedup ? it : null}
    container "bskubi/pairtools:1.0.4"

    input:
    tuple val(id), path(infile), val(dedup_params)

    output:
    tuple val(id), path("${id}_dedup.pairs.gz")

    shell:
    dedup_params = dedup_params ? dedup_params.collect
        {
            item ->
            return [
                "--output-stats": "--output-stats ${id}_dedup.stats.txt",
                "--output-dups": "--output-dups ${id}_dedup.dups.pairs.gz",
                "--output-unmapped": "--output-unmapped ${id}_dedup.unmapped.pairs.gz",
                "--output-bytile-stats": "--output-bytile-stats ${id}_dedup.bytile_stats.pairs.gz"
            ].get(item, item)
        }.join(" ") : ""
    
    cmd = "pairtools dedup --output ${id}_dedup.pairs.gz ${dedup_params} ${infile}"
    cmd
}

workflow Deduplicate {
    take:
        samples
    
    main:
        
        samples | filter{it.deduplicate} | set {deduplicate}

        JoinProcessResults(
            PairtoolsDedup,
            [deduplicate, samples],
            ["id", "latest", "dedup_params"],
            ["id", "dedup_pairs"],
            ["id"],
            false,
            "dedup_pairs"
        ) | set {samples}

        if (params.general.get("last_step") == "Deduplicate") {
            channel.empty() | set{samples}
        }

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
                condition, bioreps ->
                condition += bioreps[0].findAll{entry -> bioreps.every {map -> map[entry.key] == entry.value}}
                condition.latest = bioreps.collect{it.latest}
                condition.id = "${condition.condition}"
                condition.is_condition = true
                condition
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

        if ("BiorepsToConditions" in params.general.get("qc_after")) {
            QCReads(samples, "BiorepsToConditions")
        }

    emit:
        samples
}

process PairtoolsSelect {
    publishDir params.general.publish.select ? params.general.publish.select : "results",
               saveAs: {params.general.publish.select ? it : null}
    container "bskubi/pairtools:1.0.4"

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
    min_distances += condition.min_dist_rr != null ? ["--":condition.min_dist_rf] : [:]
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
}

workflow Select {
    take:
        samples
    
    main:        
        JoinProcessResults(
            PairtoolsSelect,
            [samples],
            ["id", "latest", "select_params", "select_condition"],
            ["id", "select_pairs"],
            ["id"],
            false,
            "select_pairs"
        ) | set{samples}

        
        if ("Select" in params.general.get("qc_after")) {
            QCReads(samples, "Select")
        }

        if (params.general.get("last_step") == "Select") {
            channel.empty() | set{samples}
        }




    emit:
        samples
}


workflow Downsample {
    take:
        samples
    
    main:
        /* Not currently in use
        */
        samples
            | map{tuple(it.subMap("is_biorep", "is_condition"), it.id)}
            | groupTuple()
            | view
    
    emit:
        samples
}

process CoolerZoomify {
    publishDir params.general.publish.mcool ? params.general.publish.mcool : "results",
               saveAs: {params.general.publish.mcool ? it : null}
    conda "cooler"
    maxForks 2

    input:
    tuple val(id), path(infile), path(chromsizes), val(pairs_format), val(assembly), val(matrix), val(make_cool), val(make_mcool)

    output:
    tuple val(id), path("${id}.mcool")

    shell:
    min_bin = matrix.resolutions.min()
    bins = matrix.resolutions.join(",")

    cmd = ["cooler cload pairs"] + make_cool +
           ["--assembly ${assembly}",
           "--chrom1 ${pairs_format.chrom1}",
           "--pos1 ${pairs_format.pos1}",
           "--chrom2 ${pairs_format.chrom2}",
           "--pos2 ${pairs_format.pos2}",
           "${chromsizes}:${min_bin}",
           "${infile} ${id}.cool",
           "&& cooler zoomify"] + make_mcool +
           ["--resolutions '${bins}'",
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
            ["id", "latest", "chromsizes", "pairs_format", "assembly", "matrix", "make_cool", "make_mcool"],
            ["id", "mcool"],
            ["id"],
            false,
            "mcool"
        ) | set{samples}


    emit:
        samples
}

process JuicerToolsPre {
    /*
        Juicer Tools Pre documentation: https://github.com/aidenlab/juicer/wiki/Pre

        We use version 1.22.01 as 2.0+ versions are in development and certain
        features available in version 1 are unavailable in 2.
    */
    publishDir params.general.publish.hic ? params.general.publish.hic : "results",
               saveAs: {params.general.publish.hic ? it : null}
    container "bskubi/juicer_tools:1.22.01"
    maxForks 2

    input:
    tuple val(id), path(infile), path(chromsizes), val(pairs_format), val(matrix)

    output:
    tuple val(id), path("${id}.hic")

    shell:
    outfile = "${id}.hic"
    min_bin = matrix.resolutions.min()
    bins = matrix.resolutions ? "-r ${matrix.resolutions.join(',')}" : ""

    cmd = ["java -Xmx20g -jar /app/juicer_tools.jar pre",
            bins,
           "${infile} ${outfile} ${chromsizes}"]
    cmd.removeAll([null])
    cmd.join(" ")
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

process MustacheDiffloops{
    publishDir "results/loops"
    container "bskubi/mustache"

    input:
    tuple val(prefix), val(id1), path(mx1), val(id2), path(mx2), val(mustache_params)

    output:
    tuple val(id1), val(id2), path("${prefix}.loop1"), path("${prefix}.loop2"), path("${prefix}.diffloop1"), path("${prefix}.diffloop2")

    shell:
    cmd = ["/bin/bash -c 'source ~/.bashrc && python mustache/mustache/diff_mustache.py",
           "-f1 ${mx1} -f2 ${mx2}",
           "-o ${prefix}"] + mustache_params
    cmd = cmd.join(" ")
    cmd += "'"
    cmd = "/bin/bash -c 'source ~/.bashrc && python mustache/mustache/diff_mustache.py --help'"
    print(cmd)
    ""
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
        | view

    emit:
    samples
}

workflow {

    channel.fromPath(params.general.samples, checkIfExists: true)
        | splitCsv(header: true)
        | map{it.id = it.sample_id; it}
        | AssignParams
        | TryDownloadMissingReferences
        | MakeMissingChromsizes
        | MakeMissingIndex
        | MakeMissingDigest
        | Align
        | Parse
        | IngestPairs
        | OptionalFragtag
        | TechrepsToBioreps
        | Deduplicate
        | BiorepsToConditions
        | Select
        | MakeHic
        | MakeMcool
        //| CallLoops
}

