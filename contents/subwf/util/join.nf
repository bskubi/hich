workflow keydiff {
    take:
        left
        right
        by
        how

    main:

    left_keys = left
        | map{it.subMap(by)}
        | unique
        | map{[it, true]}

    right_keys = right
        | map{it.subMap(by)}
        | unique
        | map{[it, true]}

    missing =
        left_keys
        | join(right_keys, remainder: true)
        | branch {
            neither: it[1] &&  it[2] //key present in both
            right:   it[1] && !it[2] //key in left only
            left:   !it[1] &&  it[2] //key in right only
        }

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

    // Get the BY-values to join on and convert to an ArrayList if it's not already
    by = condition.get("by")
    by = by instanceof List ? by : [by]

    // Get the suffix for non-BY values that are found in the left and right channels.
    suffix = condition.get("suffix", "_right")

    // Determine what type of join to do, defaulting to a left join.
    how = condition.get("how", "left")

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

    left = left | map{[it.subMap(by), it]} | groupTuple
    right = right | map{[it.subMap(by), it]}

    dominant = condition.dominant ?: "left"
    assert dominant in ["left", "right"], "In sqljoin, condition.dominant must be 'left' (default) or 'right' but was ${dominant}"

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