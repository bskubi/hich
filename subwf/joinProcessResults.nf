include {keydiff; sqljoin} from './sqljoin.nf'

def getIdx(item, index) {
    /* If item is a list, get value at index, where too-large indexes are set
       to index of final item. If item is not a list, return it.
    */
    if (item instanceof List) {
        validIndex = Math.min(index, item.size() - 1)
        return item[validIndex]
    }
    return item
}

workflow JoinProcessResults {
    /*
        proc -- the process to call
        channels -- a list of channels.
            The first channel contains single-hashmap items that have as a
            subset of their keys the required process inputs.

            Outputs from the process will be rejoined to channel[0], then
            channel[0] will be joined to channel[1], channel[1] to channel[1], ...
            and the result of the final join will be emitted.
        input -- a list of keys to extract from channel[0] hashmap items and
            use as inputs to the process
        output -- a list of keys in the order expected from the process
        join_by -- either a non-List universal join key or a List of join keys
               Iterate through the channels pairwise ([ch0, ch1], [ch1, ch2]...)
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
        condition -- if non-falsey, used to subMap channel[0] and passed as a val
                  parameter to the process
        latest -- if non-falsey, sets the 'latest' parameter of the channel items
                  used as inputs to the process to the specified output
                  from the process.

            Motivation:

            Using hashmaps to keep track of process inputs and outputs has
            numerous advantages, but adds complications.
            
            First, although Nextflow can accept hashmaps as val() params,
            it won't stage files stored in them to the work directory, making
            those files inaccessible to the process.

            Second, any changes to process inputs trigger a rerun of the process.
            We want to avoid this in cases where changes are made to hashmap
            keys that are not used by the process.

            This requires extracting the specific desired elements and passing
            them to the process, but this in turn means that the process can't
            just emit the full hashmap. We therefore have to have a way of
            extracting the desired elements, passing them to the process,
            and rejoining the process outputs to the original hashmap. This
            is what the JoinProcessResults workflow accomplishes.
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