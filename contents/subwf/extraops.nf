def isExistingFile(it) {
    // Type-agnostic way to check if file exists for any file class having an exists() method.
    return it && it.metaClass.respondsTo(it, 'exists') && it.exists()
}

def skip(step) {
    
    def excluded = params.containsKey("runOnly") && !params.runOnly.split().contains(step)
    def skipped = params.containsKey("skip") && params.skip.split().contains(step)
    return excluded || skipped
}

workflow keydiff {
    take:
        left
        right
        by
        how

    main:
    /*
        CONTEXT:

        Facilitate sql-like joins on Nextflow channels containing LinkedHashMaps

        PURPOSE:

        During a left or right join, identify missing keys that should be replaced
        by default values.

        IMPLEMENTATION:

        Take two channels LEFT and RIGHT joined on a subset of keys BY
            Channel items should be HashMaps
        Return:
            LEFT KEYS not in RIGHT
            RIGHT KEYS not in LEFT
    */
    
    /* 
        PURPOSE:
        
        Extract unique BY values from LEFT and RIGHT channels
        
        IMPLEMENTATION:
        
        Add a 'true' value as this will help keep track of whether a BY value
        was left-only, right-only, or in both.
    */
    left_keys = left
        | map{it.subMap(by)}
        | unique
        | map{[it, true]}

    right_keys = right
        | map{it.subMap(by)}
        | unique
        | map{[it, true]}

    /*
        PURPOSE:

        Split which BY-values are missing in the left (i.e. right-only), missing
        in the right (i.e. left-only), or missing in neither (i.e. found in both)

        IMPLEMENTATION:
        unique operator above means left_keys and right_keys are non-redundant

        remainder: true emits [left, null] or [null, right] when a BY value is
        left-only or right-only.

        join will result in a 3-element tuple [BY-value, inLeft, inRight]
        where inLeft and inRight are either true or null (null is the default
        Nextflow substitutes for joins when 'remainder' is true).
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
        CONTEXT:

        For left joins, return any right-only values
        For right joins, return any left-only values

        PURPOSE:

        Return the missing BY-values that are needed depending on the join type.

        IMPLEMENTATION:
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

    /*
        MOTIVATION: HashMap-based traceability

        A traceable workflow tracks every file and attribute produced for each
        sample in a single per-sample HashMap under meaningful keys. For example,
        the initial fastq files are under sample.fastq1 and sample.fastq2. When
        alignment producs a bam file, the file handle is stored at sample.sambam. 

        CONTEXT:

        Nextflow poses a few challenges for traceable workflows.
            If the entire sample HashMap is passed in as a process input,
            Nextflow will yield a different task cache for the process any time
            irrelevant Hashmap attributes change, breaking the -resume feature.

            In order to be staged, Nextflow requires paths to be extracted from
            the HashMap and passed in the format path(pathVariable) in the input
            and output sections.

        Hich solves this by assigning each sample or resource file a unique ID.
        This, along with the necessary process inputs, are extracted from the
        HashMap, then the outputs (inluding the ID) are built into a new HashMap,
        and the ID can be used as a join key to associate the process outputs
        with the original per-sample HashMap.

        The problem is that Nextflow's standard 'join' operator does not really
        implement a true SQL inner join (contrary to the documentation's claims)
        and does not work directly on HashMap keys.

        PURPOSE:

        Implement a true SQL-like inner, full, left, and right join.

        IMPLEMENTATION:

        ARGUMENTS:
        left
        right
        condition
            by
            suffix
            how
    */

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

    /*
        SETUP: Extract and format by, suffix, and how
    */

    // Get the BY-values to join on and convert to an ArrayList if it's not already
    by = condition.get("by")
    by = by instanceof List ? by : [by]

    // Get the suffix for non-BY values that are found in the left and right channels.
    suffix = condition.get("suffix", "_right")

    // Determine what type of join to do, defaulting to a left join.
    how = condition.get("how", "left")

    /*
        1. Depending on join type, determine missing keysets that would
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

    /*
        2.

        PURPOSE:

        Extract BY values in preparation for cross, grouping the LEFT hashmaps
        that have the same BY value since this is a requirement of
        Nextflow's cross operator.

        IMPLEMENTATION:

        This leaves the channels formatted as follows:

        LEFT:
        channel([BY1, [HASHMAP1, HASHMAP2, ...]], [BY2, [HASHMAP3, HASHMAP4, ...]])

        RIGHT:
        channel([BY1, HASHMAP1], [BY2, HASHMAP2])
    */

    left = left | map{[it.subMap(by), it]} | groupTuple
    right = right | map{[it.subMap(by), it]}

    dominant = condition.dominant ?: "left"
    assert dominant in ["left", "right"], "In sqljoin, condition.dominant must be 'left' (default) or 'right' but was ${dominant}"
    /*
        CONTEXT:

        At this point, LEFT items with a matching BY value are grouped as a
        single LEFT channel item (see above). We need to associate all samples
        with matching BY values and reshape them into single, joined HashMap
        items.

        IMPLEMENTATION:
        
        Nextflow's cross operator will associate each LEFT item with any RIGHT
        items having the same BY value. For example, if there are two LEFT and
        two RIGHT samples having the BY value BY1, cross yields:

            channel(
                [[BY1, [LEFT_MAP1, LEFT_MAP2]], [BY1, RIGHT_MAP1]],
                [[BY1, [LEFT_MAP1, LEFT_MAP2]], [BY1, RIGHT_MAP2]]
            )
        
        For each item, we loosely want to emit items in the format:

            [LEFT_MAP1 + RIGHT_MAP1]
            [LEFT_MAP2 + RIGHT_MAP1]
            [LEFT_MAP1 + RIGHT_MAP2]
            [LEFT_MAP2 + RIGHT_MAP2]

        The caveat is that we have to deal with a situation where a LEFT_MAP
        and its RIGHT_MAP share non-BY keys. To handle this, we declare one of
        the channels "dominant" (left by default) and a suffix to add one or
        more times to the non-dominant channel. If the suffix is blank, then
        the non-dominant channels' non-BY items are replaced by the dominant
        channel's non-BY items at the same key. Otherwise, the suffix is added
        to the non-dominant non-BY keys until they don't match any keys in the
        dominant-channel HashMap.
    */
    left
    | cross(right)
    | map{
        joinedMapList = []
        leftMapList = it[0][1]
        rightMap = it[1][1]

        // Iterate through each leftrow and combine with rightMap
        leftMapList.each() {
            leftMap ->

            joinedMap = dominant == "left" ? leftMap : rightMap

            // Iterate through each key in rightMap,
            // appending suffix if needed to avoid clashes with leftrow
            // and adding it to the joinedMap.
            (dominant == "left" ? rightMap : leftMap).each() {
                key, value ->

                // Add non-keyset keys from right to left.
                if (!by.contains(key)) {
                    
                    // Add as many suffixes as necessary to avoid clash
                    while (suffix != "" && joinedMap.containsKey(key)) {
                        key = key + suffix
                    }

                    /* Add the right key value if either criterion holds:
                        - There is no name clash, or
                        - The left table is not declared dominant
                    */
                    if (!joinedMap.containsKey(key)) {
                        joinedMap += [(key):value]
                    }
                }
            }
            // Add the new hashmap to the list of hashmaps
            joinedMapList += [joinedMap]
        }

        // List of per-keyset joins
        joinedMapList
    }
    | collect       // Currently, channel items are lists of per-keyset joins.
    | flatMap       // Collects to a single list of lists, then flattens
    | set{joined}   // Yielding the desired result of single-hashmap join items.
    
    emit:
    joined
}

def pack(addTo, addMe, by = "id") {
    sqljoin(addTo, addMe, [by:by, suffix: "", dominant:"right"])
}

def label(map, lbl) {
    return map?[(lbl)]?.toString()?.length() > 0
}

def isTechrep(map) {return label(map, "techrep") && label(map, "biorep") && label(map, "condition")}
def isBiorep(map) {return (!label(map, "techrep")) && label(map, "biorep") && label(map, "condition")}
def isCondition(map) {return (!label(map, "techrep")) && (!label(map, "biorep")) && label(map, "condition")}
def isSingleCell(map) {return map.cellBarcodeField || map.isSingleCell}

def emptyOnLastStep(step, samples) {
    def isExplicitLastStep = (params.containsKey("lastStep") && params.get("lastStep") == step)
    def isLastStep = (step == "End") || isExplicitLastStep
    def hasViewLastStep = params.containsKey("viewLastStep") && params.get("viewLastStep")
    if (isLastStep && hasViewLastStep) {
        samples
            | map {
                sample ->
                
                params.viewLastStep instanceof Boolean ? sample : sample.subMap(params.viewLastStep.split())}
            | view
    }
    return isExplicitLastStep ? channel.empty() : samples
}

// New extraops

def groupHashMap(ch, by, sortBy = ["id"]) {
    ch
    | map{tuple(it.subMap(by), it)}
    | groupTuple
    | map{it[1]}
    | map{
        mapList ->              // Sort the hashmaps to ensure same output for same set of samples on repeated runs
        if (sortBy) {
                mapList.sort {
                map1, map2 ->
                sortBy.collect {    // sortBy is a list of keys to sort by, in descending order of priority
                    key ->
                    map1[key] <=> map2[key]     // -1, 0, 1 depending on comparison outcome
                }.findResult {it != 0 ? it : null} ?: 0 // We collect all the comparisons for each sortBy key, then take the first nonzero as the sort order
            }
        } else {
            mapList
        }
    }
}

// Called from within map function on transposedSamples
def coalesce (transposedSample, defaultWhenDifferent = "_unchanged", whenDifferent = [:]) {
    // effects of whenDifferent values:
    // _shrink -> take unique values in original order of first apperance
    // _drop -> drop the key:value pair
    // _nullify -> keep the key, set value to null
    // _error -> throw an exception
    // anything else -> replace the value
    // empty values are retained unchanged

    
    def result = [:]
    transposedSample.each {
        key, value ->

        // LinkedHashSet preserves the original order of the ArrayList
        def valueSet = value instanceof ArrayList ? value as LinkedHashSet : [value] as LinkedHashSet

        // Coalesce if all entries are identical. [k:[1, 1, 1]] coalesces to [k:1]
        if (valueSet.size() == 1) {
            result += [(key):valueSet.first()]
        } else if (valueSet.size() > 1) {
            // Determine strategy for this key when value is non-homogeneous
            def todo = whenDifferent.get(key, defaultWhenDifferent)

            // Error if specified
            assert todo != 'error', "In coalesce, ${key} is set to error when values are non unique. Values are ${value} for sample:\n${sample}"

            // Otherwise, use the 'todo' strategy to decide how to alter the value
            switch (todo) {
                case "_unchanged":
                    result += [(key): value]        // Leave unaltered
                    break
                case "_shrink":
                    result += [(key): valueSet]     // Keep first instance of each distinct item in value in original order of first appearance
                    break
                case "_drop":                       // Drop it from the result
                    break
                default:
                    result += [(key): todo]         // Replace with value of 'todo'
            }

        } else if (valueSet.size() == 0) {
            result += [(key): value]                // If there is nothing in value, keep the empty list as the value
        }
    }
    return result
}

def columns (mapList, options = [:]) {
    // Get set of all keys from all maps
    def allKeys = mapList.collectMany{it.keySet()}.toSet()

    // Extract parameters to list
    def transposed = [:]

    // Iterate through all maps in the mapList
    mapList.each {
        map ->

        // Iterate through all keys from all maps
        allKeys.each {
            key ->

            /* For each key, get the current map's value or a default if supplied.
               Verify the value is non-null/missing or is OK to be null.
               Add it to the previous value list and update the transposed map.
            */
            def value = map.get(key, options.defaults?[(key)])
            if (value != null || !options.dropNull) {
                def previous = transposed.get(key, [])
                assert value != null || options.nullOK?.contains(key), "In call to 'columns' without dropNull, '${key}' is not in nullOK, and is missing/null for:\n${map}"
                def valueList = previous + [value]
                def updatedItem = [(key): valueList]
                transposed += updatedItem
            }
        }
    }

    // Return the transposed map
    transposed
}

def rows (columnsMap, options = [:]) {
    // Get the keys and the transposed values
    def keys = columnsMap.keySet().toList()
    // Nextflow has an annoying behavior where if you pass to a Path process output a list with a single item,
    // it silently converts it to a file object instead of keeping it a list of file objects. To avoid weird
    // bugs when this happens, convert non-ArrayList values to ArrayList.
    
    def values = columnsMap.values().toList()
    assert values.every {it instanceof ArrayList} || values.every {(it instanceof ArrayList && it.size() == 1) || !(it instanceof ArrayList)}
    def formattedValues = values.collect({it instanceof ArrayList ? it : [it]})
    def transposed = formattedValues.transpose()

    // Create N hashmaps, each containing the values at index i for the corresponding keys
    def result = transposed.collect { row ->
        [keys, row].transpose().collectEntries()  // Build a new map for each row
    }
    result
}

def constructIdentifier(map) {return map.subMap("condition", "biorep", "techrep", "aggregateProfileName").values().join("_")}


def configParamToHashSet(val) {
    if (val instanceof ArrayList) return new HashSet(val)
    if (val instanceof HashSet) return val
    return new HashSet([val])
}


def createCompositeStrategy(strategyKeys, strategyMap, combineHow = [:]) {

    def compositeStrategy = [:]
    subStrategies = strategyKeys ? strategyMap.subMap(strategyKeys) : [:]
    subStrategies.each {
        _, subStrategy ->
        subStrategy.each {
            key, val ->
            
            def newVals = configParamToHashSet(val)
            def oldVals = compositeStrategy.get(key, new HashSet())
            def updated = combineHow.get(key, {a, b -> b})(oldVals, newVals)
            compositeStrategy += [(key): updated]
        }
    }

    return compositeStrategy
}

def filterSamplesByStrategy(samples, strategy) {
    if (!strategy) return samples

    def reservedKeywords = ["same", "different"]
    def sampleAttributeFilter = strategy.findAll {key, value -> !(key in reservedKeywords)}
    return samples | filter {
        sample ->

        def filterOn = sample.subMap(sampleAttributeFilter.keySet())
        def passesFilter = sampleAttributeFilter.every {key, select ->
            sample.get(key) in select
        }

        passesFilter
    }
}

def combineSamplesByStrategy(samples, strategy) {
    def same = strategy.get("same", [])
    def different = strategy.get("different", [])
    def sameAndDifferent = same.intersect(different)
    if (!sameAndDifferent.isEmpty()) {
        System.err.println("Warning: In filterSamplesByStrategy, comparisons on ${sameAndDifferent} are required to be same and different, so no result is obtained")
        return channel.empty()
    }

    def filtered = filterSamplesByStrategy(samples, strategy)

    def combined = filtered
    | combine(filtered)
    | filter{it[0].id <= it[1].id}
    | unique
    | filter {
        s1, s2 ->

        def sameOK = same.every {key -> s1.get(key) == s2.get(key)}
        def differentOK = different.every {key -> s1.get(key) != s2.get(key)}
        sameOK && differentOK
 
    }
    return combined
}

def aggregateLevelLabel(sample) {
    if (isTechrep(sample)) return "techrep"
    if (isBiorep(sample)) return "biorep"
    if (isCondition(sample)) return "condition"
    return "unknown"
}

def datatypeFromExtension(path) {
    extensions = [".fastq": "fastq",
                  ".fq": "fastq",
                  ".sam": "sambam",
                  ".bam": "sambam",
                  ".pairs": "pairs",
                  ".mcool": "mcool",
                  ".hic": "hic"]
    pathString = path.toString()
    foundExtension = extensions.keySet().find {
        ext ->
        pathString.endsWith(ext) || pathString.contains("${ext}.")
    }
    return foundExtension ? extensions[foundExtension] : null
}

def parsePattern(String str, String parsePattern) {
    def patternPlaceholders = []

    def pattern = parsePattern.replaceAll(/\{([^{}]*)\}/) { match ->
        if (match[1].trim()) {
            patternPlaceholders << match[1]  // Track the placeholder name
            "(?<${match[1]}>.+?)"
        } else {
            ""  // Ignore empty placeholders
        }
    }
    
    def matcher = str =~ pattern
    def result = [:]
    
    if (matcher) {
        patternPlaceholders.eachWithIndex { placeholder, index ->
            result[placeholder] = matcher.group(index + 1)  // Retrieve group by index
        }
    }
    
    
    return result ?: null
}
