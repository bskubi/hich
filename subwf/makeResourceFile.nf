include {sqljoin} from './extraops.nf'
include {transact; pack; transpack} from './extraops.nf'
workflow MakeResourceFile {
    /*
        Motivation:
        
        There are several types of common resource files we may wish to download
        or produce from another resource file. Examples include reference genomes,
        chromsizes files, and restriction digest .bed files. Often multiple
        samples will use common reference files. This workflow identifies the
        set of common references needed by all the samples, uses the given process
        to obtain them, and then adds the file to all the samples that need it
        to avoid redundant downloading and processing.

        If the common resource files are missing, they are created. If they
        exist no special processing is done but it is inforced that they are
        of type 'file' rather than type 'string'. If no_change is required they
        are left as-is.

        exists -- channel of resource files specs that exist and where the values
        under 'file_key' should be converted to type 'file'
        missing -- channel of resource file specs that do not exist and need to be made
        no_change -- channel of resource files that should not be altered
        file_key -- list of keys where the subMap of hashmaps in 'exists' with these
        keys will be converted to type 'file'
        proc -- the Nextflow process to run to make missing files
        input -- the list of inputs expected by the process
        output -- the names of the outputs produced by the process
        join_by -- keys to join the 'made' outputs to the 'missing' items
    */
    take:
        exists
        missing
        no_change
        file_key
        proc
        input
        output
        join_by

    main:
        missing
        | map {
            // Extract input values from hashmap in order required by process
            elements = it.subMap(input).values().toList()
            tuple(*elements)
        }
        | unique    // Get unique inputs to avoid redundant processing
        | proc      // Call process
        | map{
            // Use given output names and process outputs as key:value pairs
            // in a hashmap
            proc_outputs ->
            result = [:]
            [output, proc_outputs].transpose().each {k, v -> result[k] = v}
            result
        }
        | set{made}
        
        // Join the process outputs
        missing = sqljoin(missing, made, [by: join_by, suffix: ""])

        // Make files specified by 'file_key' into files in the 'exists' channel
        // then concatenate with the made and unchanged resource files.
        result = exists
            | map{
                
                it[file_key] = file(it[file_key]);
                it
            }
            | concat(no_change)
            | concat(missing)
        
    emit:
        result
}
