include {keydiff; sqljoin} from './sqljoin.nf'

def getIdx(item, index) {
    /* If item is a list, get value at index, where too-large indexes are set
       to index of final item. If item is not a list, return it.
    */
    if (item instanceof List) {
        def validIndex = Math.min(index, item.size() - 1)
        return item[validIndex]
    }
    return item
}

workflow JoinProcessResults {
    /*
            Explanation:
            
            In Hich, think of each sample as a bundle of data and params.
            Pipeline processes add new data to the bundle. We have to decide
            what bundles to pass to which processes, pull just the needed data
            and params out of the bundles, and then put the process outputs
            back into the appropriate bundle. In Hich, bundles can be joined
            across channels using a key of choice, such as "id" or "sample_id".
            
            This is what JoinProcessResults accomplishes.

            The "bundles" in Hich are hashmaps. The workflow calling
            JoinProcessResults has already prefiltered the "bundles" and
            potentially extracted a submap (or a submap of a submap...) from
            the per-sample hashmaps.

            In JoinProcessResults, the per-sample hashmaps are submapped by
            "input" keys and submitted to the given "proc". The process outputs
            are assigned to "output" keys in a new hashmap. Typically, processes
            accept an "id" as an input and output the same "id", allowing the
            outputs to be rejoined to the original hashmap bearing the same id.

            Since there may have been multiple submapping steps to yield the
            hashmap from which process inputs are extracted, JoinProcessResults
            implements a sort of iterative join. The immediate process outputs
            are joined to the first channel in "channels", and the result is
            then joined to the second channel in "channels", and so on. The
            "join_by" param specifies the key(s) to use for each join.

            It is often useful to know the data files produced by the most
            recent process call. If "latest" is specified, then the process
            output it refers to is stored under the "latest" key in the output
            hashmap.

            JoinProcessResults allows the user to optionally extract additional
            params as process inputs via the "condition" input. If non-falsey,
            it is used to additionally submap the input hashmap and these
            values are submitted as the final inputs to the process.
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