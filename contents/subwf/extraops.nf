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

    by = condition.get("by")
    by = by instanceof List ? by : [by]
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

def hashmapdiff(ch1, ch2, by, how = "left", suffix = "__joindiff__") {
    def joined = sqljoin(ch1, ch2, [by: by, how:how, suffix:"__joindiff__"])
    def keyIsInBy = {k -> by instanceof List ? k in by : k == by}

    return joined.collect {
        hashmap ->

        def diff = [:]
        hashmap.each {
            k1, v1 ->
            if (!k1.endsWith("__joindiff__") && !keyIsInBy(k1)) {
                def k2 = k1 + "__joindiff__"
                if (!k2 in hashmap.keySet()) {
                    diff[k1] = [v1]
                } else if (hashmap[k2] != v1) {
                    v2 = hashmap[k2]
                    if (v1 instanceof Path && v2 instanceof Path) {
                        if (v1.getFileName().toString() != v2.getFileName().toString()) {
                            diff[k1] = [v1, hashmap[k2]]
                        }
                    } else {
                        diff[k1] = [v1, hashmap[k2]]
                    }
                    
                }
            }
        }
        diff
    }
}

def toHashMap(keys, vals) {
    if (!(vals instanceof ArrayList) || vals.size() == 1) {
        return [keys:vals]
    }
    return [keys, vals].transpose()
                       .collectEntries { [it[0], it[1]] }
}

def transact (proc, ch, input, output, tags, options) {
    /* Extract process inputs from hashmap channel items, call the process,
    and rebuild a hashmap with process outputs. "Tag" some outputs by assigning
    them to a second key.
    */
    
    extracted = ch
        | map {
            hashmap ->
            def proc_inputs = input.collect{
                key ->
                
                err = [
                    "In ${proc} with tags ${tags} and options ${options}, ${key} is required but is '${hashmap.get(key)}' for:",
                    "${hashmap}",
                    "\nOne possible cause is a mismatched resource file or processing parameters for a biological replicate or condition produced by a merge.",
                    "If so, you can fix this either by making all resource files and processing parameters identical for all input samples for the condition or",
                    "by specifying specific resource files and processing parameters for the biological replicate or condition in nextflow.config"
                ].join("\n")
                nullOk = options.get("nullOk")
                            ? hashmap.get(key) in options.get("nullOk") || key == options.get("nullOk")
                            : false
                assert nullOk || hashmap.get(key), err
                hashmap.get(key)
            }

            proc_inputs.size() == 1 ? proc_inputs[0] : tuple(*proc_inputs)
        }
    
    extracted = options.get("filter_unique") ? extracted | unique : extracted

    return extracted
        | proc
        | map{
            proc_outputs ->
            def result_map = toHashMap(output, proc_outputs)
            
            options.get("keep") ? result_map = result_map.subMap(options.get("keep")) : null
            options.get("remove") ? result_map -= options.get("remove") : null

            tags.each {
                tag, orig ->
                result_map[tag] = result_map[orig]
            }
            
            result_map
        }
}

def pack(channels, joinBy = "id") {
    // Extract first and remaining channels
    def first = channels[0]
    def rest = channels[1..-1]

    

    // Iteratively join the channels from first to last.
    // If joinBy is a list, join on the joinBy element corresponding to each
    // item. Otherwise, join on joinBy every time.
    sizeMatch = !(joinBy instanceof List)
                || joinBy.size() == channels.size() - 1
    assert sizeMatch, "In extraops.nf, there must be a single non-list joinBy or one joinBy for every join"
    result = rest.inject(
        first,
        {
            addMe, addTo ->
            // Get the key to join on
            idxOfAddMe = channels.indexOf(addMe)
            def by = joinBy instanceof List ? joinBy[idxOfAddMe] : joinBy            

            // Join the previous results (addMe) to the new (presumably larger)
            // hashmaps in addTo. On overlapping columns for a given
            // channel item, addMe's values replace addTo's values.
            sqljoin(addTo, addMe, [by: by, suffix: ""])
        }
    )
    return result
}

def transpack (proc, channels, input, output, tags = [:], by = "id", options = [:]) {
    // Convenience function to call transact followed by pack.
    def channelsList = channels instanceof List ? channels : [channels]
    def obtained = transact(proc, channelsList[0], input, output, tags, options)
    return pack([obtained] + channelsList, by)
}

def source (produceProc, ch, key, input, output, namer, by, needsTest) {
    // Identify channel items that need the resource file
    def needs = ch
        | branch {
            hashmap ->

            yes: needsTest(hashmap)
            no: true
        }
    
    // Give a default filename to channel items that need the resource but
    // do not have a filename already. Then filter for those having no
    // existing filename.
    def missing = needs.yes
        | map {
            hashmap ->
            
            def filename = hashmap.get(key)
            hashmap += filename ? [(key):file(filename)] : [(key):namer(hashmap)]
            def err = "During source, a hashmap needed a filename but filename given by namer function was '${filename}' at '${key}' for:\n${hashmap}"
            assert file(hashmap.get(key))?.name.length() > 0, err
            hashmap
        }
        | branch {
            hashmap ->

            yes: !file(hashmap.get(key)).exists() || file(hashmap.get(key)).getClass() == nextflow.file.http.XPath
            no: true
        }

    // Create the needed resource file and pack it into all the hashmaps that
    // were missing it.
    def made = transpack(produceProc, missing.yes, input, output, tags = [:], by = by, [filter_unique: true])

    // Concatenate the channels containing the made resources, those that
    // needed it but were not missing, and those that did not need it.
    return made | concat(missing.no) | concat(needs.no)
}

def sourcePrefix (produceProc, ch, dir, prefix, input, output, namer, by, needsTest, options = [:]) {
    // Identify channel items that need the resource file
    def needs = ch
        | branch {
            hashmap ->

            yes: needsTest(hashmap)
            no: true
        }
    
    // Give a default filename to channel items that need the resource but
    // do not have a filename already. Then filter for those having no
    // existing filename.
    def missing = needs.yes
        | map {
            hashmap ->
            
            hashmap += hashmap.get(dir) ? [(dir):file(hashmap.get(dir))] : [:]
            hashmap += hashmap.get(prefix) ? [:] : namer(hashmap)
            def err = "During sourcePrefix, ${hashmap.subMap(dir, prefix)} in:\n${hashmap}"
            assert hashmap.get(prefix), err
            hashmap
        }
        | branch {
            hashmap ->

            yes: !hashmap.get(dir) || !file(hashmap.get(dir)).exists()
            no: true
        }

    // Create the needed resource file and pack it into all the hashmaps that
    // were missing it.
    def made = transpack(produceProc, missing.yes, input, output, tags = [:], by = by, options + [filter_unique: true])
    
    // Concatenate the channels containing the made resources, those that
    // needed it but were not missing, and those that did not need it.
    return made | concat(missing.no) | concat(needs.no)
}

/*
    We can't extract the requested information from the samples and process in groovy.
    Instead we will need to do something like:
    For each profile
        Filter for samples having that profile
        subMap the samples
        collect the samples
        coalesce the samples
        add the args to the hashmap
        submap it on the complete set of args
    The caller will call the process and collect the results
*/

def parameterize(processName, samples, parameterizations, sampleKeys, inputOrder) {
    def argSets = parameterizations.get(processName)

    def resultChannels = channel.empty()

    argSets.each {
        argsName, args ->
        
        def newChannel = samples
            | filter{
                sample ->

                def hasKeys = sampleKeys.every{sample.get(it)}
                def hasProcess = sample.get(processName)
                def hasParameterization = argsName in sample.get(processName)

                if (hasProcess && !hasKeys) {
                    error [
                        "For ${processName} ${argsName}, the following sample does not have all required keys ${sampleKeys}:",
                        "${sample}",
                        "It only has ${sample.subMap(sampleKeys)}"
                    ].join("\n")
                }

                hasKeys && hasProcess && hasParameterization
            }
            | map{sample -> sample.subMap(sampleKeys)}
            | collect
            | map {
                sampleList ->

                // Coalesce the samples
                def result = [:]
                sampleList.each {
                    sample ->

                    sample.each {
                        k, v ->

                        result.get(k, []) << v
                    }
                }
                result
            }
            | map {
                sampleKeySubMap ->

                (sampleKeySubMap + args).subMap(inputOrder)
            }
        resultChannels = resultChannels.concat(newChannel)
    }

    return resultChannels
}

def label(hashmap, lbl) {
    return (hashmap.containsKey(lbl) &&
            hashmap.get(lbl) != null &&
            hashmap.get(lbl).toString().trim().length() >= 1)
}

def isTechrep(map) {return label(map, "techrep") && label(map, "biorep") && label(map, "condition")}
def isBiorep(map) {return (!label(map, "techrep")) && label(map, "biorep") && label(map, "condition")}
def isCondition(map) {return (!label(map, "techrep")) && (!label(map, "biorep")) && label(map, "condition")}

def combinations(samples, comboKeys, renameComboKeys, constKeys) {
    toCombine = samples
        | map{sample -> sample.subMap(comboKeys).values().toList().transpose().toSet()}
        | flatten
        | collate(comboKeys.size())
    
    constPart = samples
        | map{sample -> sample.subMap(constKeys)}
    
    return toCombine
        | combine(toCombine)
        | filter{it[0] > it[comboKeys.size()]}
        | map {
            sample ->

            [renameComboKeys, sample].transpose().collectEntries{[(it[0]): it[1]]}
        }
        | combine(constPart)
        | map{it[0] + it[1]}
}