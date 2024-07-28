include {MakeMissingChromsizes} from './chromsizes.nf'
include {TryDownloadMissingReferences} from './genomeReferences.nf'
include {MakeMissingIndex} from './alignerIndex.nf'
include {MakeMissingDigest} from './makeDigest.nf'

workflow AssignParams {
    take:
        samples
    
    main:
        
        
        samples
            | map {
                sample ->

                params.defaults.each {
                    k, v ->

                    // Add default params to sample hashmaps if not present
                    !(k in sample) ? sample += [(k):v] : null
                }
                params.each {
                    k, bundle ->

                    // Overwrite previous params with id-specific params
                    // Priority is to later-specified params in nextflow.config
                    sample_id = sample.id.toString()
                    bundle_ids = bundle.ids.collect{it.toString()}

                    if (bundle.containsKey("ids") && sample_id in bundle_ids) {
                        update = bundle.clone()
                        update.remove("ids")
                        sample += update
                    }
                }
                sample
            }
            | TryDownloadMissingReferences
            | MakeMissingChromsizes
            | MakeMissingIndex
            | MakeMissingDigest
            | set{samples}


    emit:
        samples
}
