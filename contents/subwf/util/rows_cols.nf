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

def rowHashmapToColumnChannel (hashMap, keyCol, valueCol) {
    // Treat key-value pairs of hashMap as row entries in a two-column table and emit in column format
    (
        channel.of(params.aggregationPlans)
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
