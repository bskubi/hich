import groovy.json.JsonOutput

def prepForJson(obj) {
    def updated = obj
    if (obj instanceof Map) {
        updated = [:]
        obj.each {
            k, v ->
            kStr = k.toString()
            updated[(kStr)] = prepForJson(v)
        }
    } else if (obj.getClass() in [List, ArrayList]) {
        updated = []
        obj.eachWithIndex {it, idx -> updated[idx] = prepForJson(it) }
    } else {
        updated = obj.toString()
    }
    return updated
}


def withLog(cmdArg, mapArg) {
    def map = mapArg + [command: cmdArg.toString()]
    // Convert to a map to avoid cyclic entries
    def preppedMap = prepForJson(map)
    def json = JsonOutput.toJson(preppedMap)
    return "${cmdArg} && cat <<'EOF' > hich_output.json\n${json}\nEOF"
}

def stubLog(stubArg, cmdArg, mapArg) {
    def map = mapArg + [command: cmdArg.toString(), stub: stubArg.toString()]
    // Convert to a map to avoid cyclic entries
    def preppedMap = prepForJson(map)
    def json = JsonOutput.toJson(preppedMap)
    return "${stubArg} && cat <<'EOF' > hich_output.json\n${json}\nEOF"
}

def isExistingFile(it) {
    // Type-agnostic way to check if file exists for any file class having an exists() method.
    return it && it.metaClass.respondsTo(it, 'exists') && it.exists()
}

def skip(step) {
    /*
        Users may want to skip some steps, such as QC or forming a particular kind of contact matrix,
        or run only certain steps. This uses both params to define a list of steps to be skipped
        (the intersection of skip and runOnly's complement).
    */
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
    /*
        A big feature of Hich is the ability to store sample attributes in hashmaps.
        This requires extracting just the features specifically relevant to a particular
        process and passing them in. Otherwise, irrelevant changes to the hashmap would trigger
        reruns of the process. However, this then requires rejoining the outputs from the process
        to the original hashmap. We therefore pass in a unique id for each sample. The sqljoin
        function does most of the heavy lifting but is more general than required, so "pack"
        provides a simpler gloss over this needed functionality.

        Update the "addTo" channel of hashmaps with the "addMe" channel of hashmaps
        dominant is "right" because this makes addMe overwrite corresponding keys
        in addTo rather than the reverse. Suffix is blank so we overwrite rather than
        adding new keys when addTo and addMe share keys.
    */
    sqljoin(addTo, addMe, [by:by, suffix: "", dominant:"right"])
}

def label(map, lbl) {
    // Return whether map.lbl contains a non-empty string, used below to determine
    // if the techrep, biorep and condition keys are present and specified
    return map?[(lbl)]?.toString()?.length() > 0
}

/*
    Determine if the techrep, biorep, condition fields are uniquely specified

    For isSingleCell, the user has to specify a cellBarcodeField (which tag in a sam/bam
    file holds the cell barcode in order to extract it to the pairs cellID field)
    or has to specify isSingleCell (which can be used for pairs files where the cellID field
    is already extracted).
*/
def isTechrep(map) {return label(map, "techrep") && label(map, "biorep") && label(map, "condition")}
def isBiorep(map) {return (!label(map, "techrep")) && label(map, "biorep") && label(map, "condition")}
def isCondition(map) {return (!label(map, "techrep")) && (!label(map, "biorep")) && label(map, "condition")}
def isSingleCell(map) {return map.cellBarcodeField || map.isSingleCell}

/*
    emptyOnLastStep is used to control behavior when the workflow reaches the step specified by the --lastStep
    param, if specified. It should get passed the main samples channel and returns an empty channel if it's the 
    last step (halting execution) or the original samples channel otherwise.

    If --viewLastStep is specified it will display the contents of the samples channel after the last step finishes.
    If a space-separated string is specified, it will turn that into a list and submap just those sample attributes for viewing.
*/
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

def groupHashMap(ch, by, sortBy = ["id"]) {
    /* Essentially a groupby operation on a channel of hashmaps

        Example:

        [[param1: 1, param2: 1], [param1: 1, param2: 2], [param1: 2, param2: 3], [param1: 2, param2: 4]]

        This is grouped to:

        [param1: [1, 1], param2: [1, 2]]
        [param1: [2, 2], param2: [3, 4]]

        Arguments:
            ch -- the sample hashmap channel
            by -- the list of attributes to group by
            sortBy -- typically just use id. This function is often used to convert samples to columnar format where a single
            hashmap containing lists of sample attributes under each attribute name is provided to the process. Changes in order
            typically should not affect process behavior, but will trigger a Nextflow rerun unnecessarily. We therefore use sortBy
            to sort the samples to prevent these unnecessary reruns.
    */
    ch
    // Group maps to the format [[bySubmap], [list of samples]] and extract [list of samples]
    // Note that because the order in which samples are returned is unpredictable, this will result in a different
    // order of samples on different runs, which is why samples are sorted later.
    | map{tuple(it.subMap(by), it)} 
    | groupTuple
    | map{it[1]}

    // Sort the hashmaps if requested
    | map{
        mapList ->

        // Sort the hashmaps to ensure same output for same set of samples on repeated runs
        if (sortBy) {
                mapList.sort {
                map1, map2 ->

                // sortBy is a list of keys to sort by, in descending order of priority
                // We collect all the comparisons for each sortBy key, then take the first nonzero as the sort order
                sortBy.collect {
                    key ->

                    // -1, 0, 1 depending on comparison outcome
                    map1[key] <=> map2[key]
                }.findResult {it != 0 ? it : null} ?: 0
                
            }
        } else {
            mapList
        }
    }
}

// Called from within map function on transposedSamples
def coalesce (transposedSample, defaultWhenDifferent = "_unchanged", whenDifferent = [:]) {
    /* Combine lists of values that are identical into a single non-list value

        Example:

        [param1: [1, 1, 1], params: [2, 3, 4]] is converted to [param1: 1, params2: [2, 3, 4]]

        Arguments
            transposedSample -- The sample in columnar format
            defaultWhenDifferent: [str] -- Determines how values that are different are treated
                _shrink -> take unique values in original order of first apperance
                _drop -> drop the key:value pair
                _nullify -> keep the key, set value to null
                _error -> throw an exception
                anything else -> replace the value
                empty values are retained unchanged

            whenDifferent: map[str, str] -- Key-specific behavior when different. The values
            can be any of the appropriate values for defaultWhenDifferent and are applied to
            specific keys in transposedSample. For a particular key, whenDifferent's value will be
            used if the key is present in whenDifferent, and defaultWhenDifferent otherwise.
    */

    def result = [:]

    // Iterate through values in the transposed sample
    transposedSample.each {
        key, value ->

        // LinkedHashSet preserves the original order of the ArrayList
        def valueSet = value instanceof ArrayList ? value as LinkedHashSet : [value] as LinkedHashSet

        if (valueSet.size() == 1) {
            // If there's only one value, just use that.
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
    def allKeys = mapList.collectMany{
        it.keySet()
    }.toSet()

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

def mapToColumns (hashMap, keyCol, valueCol) {
    [
        (keyCol): hashMap.keySet().toList(), 
        (valueCol): hashMap.values().toList()
    ]
}

def rows (columnsMap) {
    /*
        Convert columnsMap from columnar to row format

        In columnar format, each parameter is a key in columnsMap, with a list
        of values corresponding at a given index to a particular sample.

        In row format, we have a row of maps where each map corresponds to a sample,
        with keys the parameter names and values the parameter values.

        Columnar example:
            [
                param1: [1, 2, 3]
                param2: ["a", "b", "c"]
            ]

        Row example:
            [
                [param1: 1, param2: "a"],
                [param1: 2, param2: "b"],
                [param1: 3, param2: "c"],
            ]

        Arguments
            columnsMap: map[str, ArrayList | Any] -- A channel containing key names
            associated with lists of per-sample values to be converted to a row format
    */

    // Get the keys and the transposed values
    def keys = columnsMap.keySet().toList()  
    def values = columnsMap.values().toList()

    // Ensure all values are lists
    def formattedValues = values.collect({it instanceof ArrayList ? it : [it]})

    // This converts the column format of the values to row format
    // https://docs.groovy-lang.org/latest/html/groovy-jdk/java/util/List.html#transpose()
    def transposed = formattedValues.transpose()

    // Create N hashmaps, each containing the values at index i for the corresponding keys
    def result = transposed.collect {
        row ->
        // Keys and row are equal-length lists which are transposed into a list of individual
        // key value pairs like [[key1: value1], [key2: value2]]. collectEntries converts
        // this to a single map like [key1: value1, key2: value2]
        [keys, row].transpose().collectEntries()
    }

    result
}

/*
    Hich depends on each sample having a unique sample.id attribute for joining process results to the
    appropriate sample hashmap. A legible name is also convenient for troubleshooting. If the user
    wants to let Hich build unique ids automatically, they should specify unique conditions, bioreps and techreps
    and not use the _ character in order to ensure that all ids will be unique. The aggregateProfileName is also
    included because new copies of the input samples are produced for each aggregateProfile.
*/
def constructIdentifier(map) {
    return map.subMap("condition", "biorep", "techrep", "aggregateProfileName").values().join("_")
}


def asHashSet(val) {
    // Convert a non-HashSet val into a HashSet
    if (val instanceof ArrayList) return new HashSet(val)
    if (val instanceof HashSet) return val
    return new HashSet([val])
}


def createCompositeStrategy(strategyKeys, strategyMap, combineHow = [:]) {
    /* A composite strategy is a hashmap in which keys are sample attributes and values are lists of permitted sample attribute values. It is created by combining one or more individual strategies specified in params.sampleSelectionStrategies.

        strategyKeys -- keys in strategyMap for the sub-strategies to combine (i.e. analysisPlan)
        strategyMap -- [strategyKey: selectionStrategy] map-of-maps, typically params.sampleSelectionStrategies
        combineHow -- not currently used, but permits defining how to combine strategies when
            there are conflicts by passing a [key: closure] map where
            closure(oldVals, newVals) outputs the updated value for the key. By
            default, the later-specified strategy has precedence.

        Returns empty map if no keys or strategies are supplied.
    */

    def compositeStrategy = [:]

    // Extract the values associated with individual selected strategies to form the composite strategy
    def subStrategies = strategyKeys ? strategyMap.subMap(strategyKeys).values() : []

    subStrategies.each {
        subStrategy ->

        subStrategy.each {
            key, val ->

            // Converts val from a single element or ArrayList into a HashSet of elements
            def newVals = asHashSet(val)

            // Handle situations where the same key is defined in more than one composite strategy
            // combineHow may contain a per-key method to define how to do the replacement.
            // The default behavior is to replace values from earlier-specified keys with
            // newly-specified keys. For example, if the keys are ["strategy1", "strategy2"]
            // and both strategies have a value "v1", then v1 will take the value for strategy2 by default.
            def oldVals = compositeStrategy.get(key, new HashSet())
            def updated = combineHow.get(key, {a, b -> b})(oldVals, newVals)

            // Add the value to the composite strategy
            compositeStrategy += [(key): updated]
        }
    }

    return compositeStrategy
}

def filterSamplesByStrategy(samples, strategy) {
    /*  After a composite strategy is built, filter for samples for which all sample attributes are present and are in the list of permitted values specified by the composite strategy.

        samples - a channel of sample hashmaps
        strategy - a hashmap as [attributeName: [permittedValues]]
    */

    // Return all samples if no strategy is specified
    if (!strategy) return samples

    // Remove reserved keywords from the set of sample-specific strategies
    def reservedKeywords = ["same", "different"]
    def sampleAttributeFilter = strategy.findAll {key, value -> !(key in reservedKeywords)}

    def filtered = samples | filter {
        sample ->

        def passesFilter = sampleAttributeFilter.every {key, select ->
            // Check that the sample attribute value is in the whitelisted values specified in the composite strategy
            sample.get(key) in select
        }

        passesFilter
    }

    // For each sample, collect into a single list and ensure that at least one sample was selected
    filtered
        | collect
        | {
            assert it.size() > 0, "Error: In filterSamplesByStrategy with sample selection strategy ${strategy}, no samples matched this filter."
        }

    return filtered
}

def groupSamplesByStrategy(samples, strategy) {
    /* Get all samples having matching values of strategy.same
    */
    return samples
        | map{tuple(it.subMap(strategy.get("same", [])), it)}
        | groupTuple
        | map{it[1]}
}

def pairSamplesByStrategy(samples, strategy) {
    /*
        The diffloops workflow has to define a way to pair up samples, especially
        by enforcing that certain attributes are the same or different.

        samples -- channel of sample hashmaps
        strategy -- a composite sample selection strategy [attribute: permittedValues]
            two possible attributes are
                same: a list of attributes which must be the same to pair two samples
                different: a list of attributes which must differ to pair two samples
            other attributes are used for filtering individual samples
    */

    // Ensure there's no conflict between "same" and "different"
    def same = strategy.get("same", [])
    def different = strategy.get("different", [])
    def sameAndDifferent = same.intersect(different)
    if (!sameAndDifferent.isEmpty()) {
        System.err.println("Warning: In filterSamplesByStrategy, comparisons on ${sameAndDifferent} are required to be same and different, so no result is obtained")
        return channel.empty()
    }

    // Filter individual samples before forming pairs of samples
    def filtered = filterSamplesByStrategy(samples, strategy)

    // Obtain pairs of samples matching the "same" and "different" criteria
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

    combined
    | collect
    | map{
        assert it.size() > 0, "In pairSamplesByStrategy with sample selection strategy ${strategy}, no samples were paired, likely due to 'same' or 'different' filter."
    }

    return combined
}

def aggregateLevelLabel(sample) {
    // Return a magic string for the aggregateLevel
    if (isTechrep(sample)) return "techrep"
    if (isBiorep(sample)) return "biorep"
    if (isCondition(sample)) return "condition"
    return "unknown"
}

def datatypeFromExtension(path) {
    /*
        Look for various known extensions to extract the datatype implicitly from the
        input file so that Hich can ingest intermediate file formats appropriately
        without explicit specification by the user. This is especially helpful in
        permitting the user to use globs at the command line to feed files into Hich.
    */
    def extensions = [".fastq": "fastq",
                  ".fq": "fastq",
                  ".sam": "sambam",
                  ".bam": "sambam",
                  ".pairs": "pairs",
                  ".mcool": "mcool",
                  ".hic": "hic"]
    def pathString = path.toString()
    def foundExtension = extensions.keySet().find {
        ext ->
        pathString.endsWith(ext) || pathString.contains("${ext}.")
    }
    return foundExtension ? extensions[foundExtension] : null
}

def parsePattern(String str, String parsePattern) {
    /*
        Used to extract sample attributes from filenames, such as condition, biorep, and techrep,
        via a syntax similar to that offered by Python's parse library. In parse, users can
        extract substrings into a map with patterns like: "{condition}_{biorep}_{techrep}.fastq",
        which would take a string like "cond1_br1_tr1.fastq" and return ["condition": "cond1", "biorep": "br1", "techrep": "tr1"].
        This is easier to specify at the command line than a regex but AFAIK has no Groovy equivalent.

        This function implements this parsing functionality, returning the extracted map.
    */

    def patternPlaceholders = []

    // This regex searches the parsePattern string (i.e. "{condition}_{biorep}_{techrep}.fastq")
    // for the placeholders between braces (i.e. condition, biorep, techrep) and adds them
    // to the list of patternPlaceholders to become keys in the output map.

    // It also yields in "pattern" the list of matchers to look for in the input string "str"
    def pattern = parsePattern.replaceAll(/\{([^{}]*)\}/) { match ->
        if (match[1].trim()) {
            patternPlaceholders << match[1]  // Track the placeholder name
            "(?<${match[1]}>.+?)"
        } else {
            ""  // Ignore empty placeholders
        }
    }
    
    // This extracts the patterns from str
    def matcher = str =~ pattern

    // Combine the patternPlaceholders with the corresponding matches from "str" into an output map "result"
    def result = [:]
    if (matcher) {
        patternPlaceholders.eachWithIndex { placeholder, index ->
            result[placeholder] = matcher.group(index + 1)  // Retrieve group by index
        }
    }
    
    return result ?: null
}

def formatArg(pattern, object, sep) {
    /*
        Some processes receive either a list of values or a single non-list element
        as parameter values, but need to call a CLI command passing a delimiter-separated
        list of the received values. Other times they get nothing and should
        not pass an argument for that parameter at all. This facilitates this interconversion
        and returns an empty string if the object passed was falsey.

        NOTE this is a potential issue if the goal is to pass a boolean false
        to the CLI command, but I don't think Hich currently does this...

        pattern -- the string pattern to format the results into, like "--numbers {commaSeparatedNums}"
        object -- the element or list of elements to join (where necessary) into a delimiter-separated list
        sep -- the delimiter, like ","
    */
    // Put non-lists into a list so when join is called, it has something to (silently) operate on
    def listed = (object instanceof List || object == null) ? object : [object]
    def joined = listed ? listed.join(sep) : listed

    // Format the string with the result or return 
    return joined ? String.format(pattern, joined) : ""
}