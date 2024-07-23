workflow AssignParams {
    take:
        samples
    
    main:
        
        samples
            | map {
                sample ->
                params.defaults.each {
                    k, v ->
                    !(k in sample) ? sample += [(k):v] : null
                }
                params.each {
                    k, bundle ->
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
            | set{samples}
    emit:
        samples
}
