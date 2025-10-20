include {update} from './update.nf'

workflow runProcess {
    take:
    process_name
    PROCESS
    samples

    main:
    samples
        | map{HichEngine.prepareProcessInputs(process_name, it)}
        | PROCESS
        | map{
            output ->
            sample = HichEngine.formatProcessOutputs(process_name, output)
            if (!sample.containsKey("id")) {
                throw new Exception("Output from '${process_name}' lacked 'id' after formatting.")
            }
            [sample.id, sample]
        }
        | set{outputs}
    
    update(samples, outputs)
        | set{samples}

    emit:
    samples
}