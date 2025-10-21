class HichEngine {
    static Map getProcess(String processKey) {
        def PROCESS_INTERFACES = HichContract.PROCESS_INTERFACES
        if (!PROCESS_INTERFACES.containsKey(processKey)) {
            throw new Exception("FATAL: HichContract.PROCESS_INTERFACES does not contain key '${processKey}'")
        }
        def processMap = HichContract.PROCESS_INTERFACES[processKey]
        if (!Map.isInstance(processMap)) {
            throw new Exception("FATAL: HichContract.PROCESS_INTERFACES[${processKey}] is not a Map.")
        }
        return processMap
    }

    static Map prepareSample(String processKey, Map sample) {
        try {
            def processMap = getProcess(processKey)
            def errors = []

            def PRE_VALIDATORS = HichContract.Schema.PRE_VALIDATORS
            if (processMap.containsKey(PRE_VALIDATORS)) {
                errors = validate(processMap[PRE_VALIDATORS], sample)
            }
            if (!errors.isEmpty()) {
                throw new Exception("Pre-validation failed on process '${processKey}': ${errors.join('\n')}\nSample: ${sample}")
            }
            def TRANSFORMERS = HichContract.Schema.TRANSFORMERS
            if (processMap.containsKey(TRANSFORMERS)) {
                sample = transform(processMap[TRANSFORMERS], sample)
            }

            def POST_VALIDATORS = HichContract.Schema.POST_VALIDATORS
            if (processMap.containsKey(POST_VALIDATORS)) {
                errors += validate(processMap[POST_VALIDATORS], sample)
            }
            if (!errors.isEmpty()) {
                throw new Exception("Post-validation failed on process '${processKey}': ${errors.join('\n')}")
            }
        } catch(Exception e) {
            def errors = "Failed to process interface '${processKey}' on sample ${sample}.\n${e}"
            throw new Exception(errors)
        }
        return sample
    }

    static Map formatProcessOutputs(String processKey, List outputs) {
        def processMap = getProcess(processKey)
        def OUTPUTS = HichContract.Schema.OUTPUTS
        
        if (!processMap.containsKey(OUTPUTS)) {
            throw new Exception("No OUTPUTS defined in HichContract for process '${processKey}'.")
        }
        else {
            def schema = processMap[OUTPUTS]
            def count = schema.size()
            def actual_count = outputs.size()
            if (count != actual_count) {
                throw new Exception("In HichContract, ${processKey}.OUTPUTS expects ${count} outputs, got ${actual_count}.")
            }
            def result = [:]
            for (int i = 0; i < count; i++) {
                def output_i = outputs[i]
                def schema_i = schema[i]
                if (!schema_i.containsKey("key")) {
                    throw new Exception("In HichContract, ${processKey}.OUTPUTS contains schema ${schema_i} which is missing 'key'.")
                }
                def sample_key = schema_i.key
                def sample_value = HichWorkflowAdapter.convertToSchema(output_i, schema_i)
                result += [(sample_key) : sample_value]
            }
            return result
        }
    }

    static List formatInputs(String processKey, Map sample) {
        def processMap = getProcess(processKey)
        def INPUTS = HichContract.Schema.INPUTS

        if (!processMap.containsKey(INPUTS)) {
            throw new Exception("No INPUTS defined in HichContract for process '${processKey}'.")
        }
        else {
            def inputs = processMap[INPUTS]
            def fromSample = inputs.findAll{!it.containsKey("fromProcess") || !it.fromProcess}
            return fromSample.collect {
                input ->
                HichWorkflowAdapter.extractNextflowProcessInput(sample, input)
            }
        }
    }

    static List prepareProcessInputs(String processKey, Map sample) {
        sample = HichEngine.prepareSample(processKey, sample)
        return formatInputs(processKey, sample)
    }

    static HichExecutionContext getContext(String process_key, List input) {
        def interfaces = HichContract.PROCESS_INTERFACES
        
        if (!interfaces) {
            throw new Error("HichContract.PROCESS_INTERFACES does not exist.")
        }
        if (!Map.isInstance(interfaces)) {
            throw new Error("HichContract.PROCESS_INTERFACES is not a Map.")
        }
        if (!interfaces.containsKey(process_key)) {
            throw new Error("HichContract.PROCESS_INTERFACES.${process_key} does not exist.")
        }
        def process_interface = interfaces[process_key]

        if (!Map.isInstance(process_interface)) {
            throw new Error("HichContract.PROCESS_INTERFACES.${process_key} is not a Map.")
        }

        def COMMANDS = HichContract.Schema.COMMANDS
        def INPUTS = HichContract.Schema.INPUTS

        if (!process_interface.containsKey(COMMANDS)) {
            throw new Error("HichContract.PROCESS_INTERFACES.${process_key}.COMMANDS does not exist.")
        }
        if (!process_interface.containsKey(INPUTS)) {
            throw new Error("HichContract.PROCESS_INTERFACES.${process_key}.INPUTS does not exist.")
        }

        def context_builder = process_interface[COMMANDS].method
        def input_schema = process_interface[INPUTS]

        def expected_input_arg_types = input_schema.collect{
            if (!it.type) {
                throw new Exception("In PROCESS_INTERFACES[${process_key}].INPUTS, required 'type' key not specified for ${it}")
            }
            it.type
        }
        def expected_arg_types = [Map] + expected_input_arg_types

        def error_message
        if (!HichExecutionContextBuilder.metaClass.getMetaMethod(context_builder, expected_arg_types as Class[])) {
            error_message = "FATAL: Mismatch between HichContract schema for PROCESS_INTERFACES[${process_key}].INPUTS and HichExecutionContextBuilder.${context_builder} definition.\n"
            
        }

        def commands = process_interface[COMMANDS].commands
        def builder_input = [commands, *input]
    
        try {
            return HichExecutionContextBuilder."$context_builder"(*builder_input)
        } catch(Exception e) {
            throw new Exception("${error_message}\n${e}")
        }
    }

    static List validate(Map validation_rules, Map sample) {
        def all_errors = []
        
        validation_rules.each {
            rule, list_of_params ->
            
            if (!List.isInstance(list_of_params)) {
                throw new Exception("FATAL: Validation rule ${rule}'s value is class ${list_of_params.class.simpleName} rather than required List")
            }
            list_of_params.each {
                params ->
                if (!Map.isInstance(params)) {
                    throw new Exception("FATAL: List of params for validation rule ${rule} contains item of type ${params.class.simpleName}. All items must be Maps.")
                }
                try {
                    all_errors += HichValidator."$rule"(sample, params)
                } catch (Exception e) {
                    def expected_arg_types = [Map, Map] as Class[]
                    def params_class = params.getClass()
                    def actual_arg_types = [Map, params_class]
                    def explanation
                    if (!HichValidator.metaClass.getMetaMethod(rule, expected_arg_types)) {
                        explanation = "FATAL: HichContract contains misnamed rule '${rule}'.\n"
                    } else if (!HichValidator.metaClass.getMetaMethod(rule, actual_arg_types)) {
                        explanation = (
                            "At rule '${rule}', contract contains " +
                            "params item '${params}' of invalid type " +
                            "${params_class} (requires Map).\n"
                        )
                    }
                    def error_message = "${explanation}${e}"
                    throw new Exception(error_message)
                }
                
            }
        }
        return all_errors
    }

    static Map transform(Map transformation_rules, Map sample) {
        def all_errors = []
        
        transformation_rules.each {
            rule, list_of_params ->
            
            if (!List.isInstance(list_of_params)) {
                throw new Exception("FATAL: Transformation rule ${rule} is class ${list_of_params.class.simpleName} rather than required List")
            }
            list_of_params.each {
                params ->
                if (!Map.isInstance(params)) {
                    throw new Exception("FATAL: List of params for transformation rule ${rule} contains item of type ${params.class.simpleName}. All items must be Maps.")
                }
                try {
                    all_errors += HichTransformer."$rule"(sample, params)
                } catch (Exception e) {
                    def expected_arg_types = [Map, Map] as Class[]
                    def params_class = params.getClass()
                    def actual_arg_types = [Map, params_class]
                    def explanation
                    if (!HichTransformer.metaClass.getMetaMethod(rule, expected_arg_types)) {
                        explanation = "FATAL: HichContract contains misnamed rule '${rule}'.\n"
                    } else if (!HichTransformer.metaClass.getMetaMethod(rule, actual_arg_types)) {
                        explanation = (
                            "At transformation rule '${rule}', contract contains " +
                            "params item '${params}' of invalid type " +
                            "${params_class} (requires Map).\n"
                        )
                    }
                    def error_message = "${explanation}${e}"
                    throw new Exception(error_message)
                }
            }
        }
        return sample
    }
}