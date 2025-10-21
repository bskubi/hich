include {update} from './update.nf'
import hich.Engine

workflow runProcess {
    take:
    process_name
    PROCESS
    samples

    main:
    samples
        | map{Engine.getInput(process_name, it)}
        | PROCESS
        | set{output_channels}
    
    output_channels.result
        | map{
            output ->
            sample = Engine.formatOutput(process_name, output)
            if (!sample.containsKey("id")) {
                throw new Exception("Output from '${process_name}' lacked 'id' after formatting.")
            }
            [sample.id, sample]
        }
        | set{outputs}
    
    update(samples, outputs)
        | set{samples}

    output_channels.execution_context
        | set{execution_context}

    emit:
    samples
    execution_context
}