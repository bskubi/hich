include {asHashSet} from "./dataStructures.nf"

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