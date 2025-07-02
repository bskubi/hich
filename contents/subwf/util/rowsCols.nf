include {checkMapTypes} from './validation.nf'

def columns (List<Map> mapList, Map options = [:]) {
    /* Convert list of maps to single map containing list of values for each key
       
       Example
           mapList [[a: 1, b: 2, c: 5, d: 6], [a: 3, b: 4]]
           options [c: 5, nullAllOK: true]
           returns [a: [1, 3], b: [2, 4], c: [5, 5], d: [6, null]]

           mapList [[a: 1], [:]]
           options [dropAllNull: true]
           returns [a: [1]]

       options:
        "defaults" (Map): Default values for missing keys (defaults to [:])
        "dropNull" (List): Keys whose values will be dropped if null (defaults to [:])
        "dropAllNull" (boolean): If true, drop all null values (defaults to false)
        "nullOK" (List): Keys whose values are allowed to be null (error otherwise, defaults to [:])
        "nullAllOK" (boolean): If true, all values may be null (defaults to true if dropAllNull is true, false otherwise)
    */

    // Input validation
    def defaults = options.get("defaults", [:])
    def dropNull = options.get("dropNull", [])
    def nullOK = options.get("nullOK", [])
    def dropAllNull = options.get("dropAllNull", false)
    def nullAllOK = options.get("nullAllOK", dropAllNull)
    def validOptions = [
        "defaults": Map, 
        "dropNull": List, 
        "nullOK": List, 
        "dropAllNull": Boolean, 
        "nullAllOK": Boolean
    ]
    checkMapTypes(options, validOptions)

    // Get set of all keys from all maps
    def allKeys = mapList.collectMany{
        it.keySet()
    }.toSet()

    // Will store column-formatted mapList 
    def result = [:]

    // Iterate through all maps in the mapList
    mapList.each {
        map ->

        // Iterate through all keys from all maps
        allKeys.each {
            key ->

            def value = null

            def hasKey = map.containsKey(key)
            def hasDefault = defaults && defaults.containsKey(key)
            def nullFailure = !nullAllOK && !nullOK.contains(key)

            if (hasKey) {
                // Use supplied value
                value = map.get(key)
            }
            else if (hasDefault) {
                // Use default for missing value
                value = defaults.get(key)
            }
            else if (nullFailure) {
                // Missing value with no default, null not OK
                error("In call to 'columns' without dropNull (options = ${options}, key = ${key}, value = ${value}), '${key}' is not in nullOK, and is missing/null for:\n${map}")
            }
            
            def keepIfNull = !dropAllNull && !dropNull.contains(key)
            def keepValue = value || keepIfNull

            if (keepValue) {
                // Add to list of values for current key
                def values = result.get(key, []) + [value]
                result += [(key): values]
            }
        }
    }

    // Return the column-format map
    result
}

def rows (columnsMap) {
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

workflow columnsToRows {
    take:
    samples

    main:

    samples
    | map{rows(it)}
    | flatten
    | set{samples}

    emit:
    samples
}

def rowHashmapToColumnChannel (hashMap, keyCol, valueCol) {
    // Treat key-value pairs of hashMap as row entries in a two-column table and emit in column format
    (
        channel.of(params.aggregate)
        | map {
            [
                (keyCol): hashMap.keySet().toList(), 
                (valueCol): hashMap.values().toList()
            ]
        }
    )
}

def rowHashmapToRowChannel (hashMap, keyCol, valueCol) {
    // Treat key-value pairs of hashMap as row entries in a two-column table and emit in row format
    rowHashmapToColumnChannel(hashMap, keyCol, valueCol) | map{rows(it)} | flatten
}

