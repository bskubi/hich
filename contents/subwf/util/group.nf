include {columns; sortMapList} from './rowsCols.nf'

workflow GroupAsColumns {
    take:
    chan
    groupBy
    sortBy
    columnsOptions

    main:
    /*
        Group map channel by keys and convert each group to column format

        Example:
            chan channel.of([a: 1, b: 1, id: 2], [a: 1, b: 2, id: 1], [a: 2, b: 3, id: 3])
            groupBy ["a"]
            sortBy ["id"]
            columnsOptions [:]

            returns
                [
                    [a: [1, 1], b: [2, 1], id: [1, 2]],
                    [a: [2], b: [3], id: [3]]
                ]
        
        List groupBy: Keys to group maps by
        List sortBy: Keys to sort within groups for deterministic group outputs
        Map columnsOptions: Map of options to 'columns' method 
    */
    chan
        | map{tuple(it.subMap(groupBy), it)}
        | groupTuple                    
        | map{it[1]}
        | map{ mapList -> sortBy ? sortMapList(mapList, sortBy) : mapList }
        | map{ mapList -> columns(mapList, columnsOptions)}
        | set{chan}

    emit:
    chan
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
            assert todo != 'error', "In coalesce, ${key} is set to error when values are non unique. Values are ${value} for sample:\n${transposedSample}"

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