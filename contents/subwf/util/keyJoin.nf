include {missingKeys} from './missingKeys.nf'
include {validateMap} from './validation.nf'

workflow keyJoin {
    take:
    left
    right
    options

    main:
    
    // 1. Input validation
    // Must contain 'by', and is limited to by, suffix, and how options
    if (options instanceof String) {
        options = [by: [options]]
    }
    else if (options instanceof List) {
        options = [by: options]
    }
       
    validOptions = [
        how: String,
        by: null,
        suffix: String
    ]
    validateMap(options.subMap("by"), validOptions.subMap("by"), ["requireKeys"])
    validateMap(options.findAll{it.key != "by"}, validOptions.subMap(["suffix", "how"]), ["limitKeys", "limitTypes"])

    
    how = options.get("how", "left")
    validHow = ["left", "right", "full"]
    assert validHow.contains(how), "In sqljoin, how must be one of ${validHow} but was ${how}"

    
    // Get the join by values and convert to list if it's not already
    by = options.get("by")
    by = by instanceof List ? by : [by]

    // Get the suffix for keys other than those to join by
    suffix = options.get("suffix", "_right")

    
    if (how == "full" || how == "left") {
        // Add any keys in left that are missing in right to right.
        missingKeys(left, right, by) | set{missingFromRight}
        right | concat(missingFromRight) | set{right}
    }
    if (how == "full" || how == "right") {
        // Add any keys in left that are missing in left to left.
        missingKeys(right, left, by) | set{missingFromLeft}
        left | concat(missingFromLeft) | set{left}
    }
    
    // Ensure that all items in left and right contain all by keys
    left
        | concat(right)
        | map{
            sample ->
            by.every {
                key ->
                assert sample.containsKey(key), "In sqljoin with options ${options}, a sample is missing the join by key '${key}':\n${sample}"
            }
        }

    //Group left and right by common keys
    left
        | map{[it.subMap(by), it]} 
        | groupTuple 
        | set{left}
    right 
        | map{[it.subMap(by), it]} 
        | set{right}

    preserved = options.preserved ?: "left"
    assert preserved in ["left", "right"], "In sqljoin, options.preserved must be 'left' (default) or 'right' but was ${preserved}"

    

    left
    | cross(right)
    | map{
        joinedMapList = []
        leftMapList = it[0][1]
        rightMap = it[1][1]

        // Iterate through each leftrow and combine with rightMap
        leftMapList.each() {
            leftMap ->

            (preservedMap, otherMap) = (preserved == "left" ? [leftMap, rightMap] : [rightMap, leftMap])
            
            // Iterate through each key in rightMap,
            // appending suffix if needed to avoid clashes with leftrow
            // and adding it to the joinedMap.
             
            otherMap.each() {
                key, value ->

                // On key collisions, try renaming the key from the other map
                if (!by.contains(key)) {
                    if (preservedMap.containsKey(key)) {
                        renamedKey = key + suffix
                        assert !otherMap.containsKey(renamedKey), "Key '${renamedKey}' renamed from '${key}' in map\n${otherMap}\nalready exists in that map"
                        assert !preservedMap.containsKey(renamedKey), "Key '${renamedKey}' renamed from '${key}' in map\n${otherMap}\nalready exists in ${preserved} map ${preservedMap}" 
                        key = renamedKey
                    }

                    // If we get here, we can add the key to the preserved map
                    preservedMap += [(key):value]
                }
            }
            // Add the new hashmap to the list of hashmaps
            joinedMapList += [preservedMap.clone()]
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