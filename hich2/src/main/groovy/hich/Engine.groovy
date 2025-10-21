package hich
import hich.specs.HichSpec
import hich.specs.HichSpec.SpecKey
import hich.specs.plans.keys.PlanKey
import hich.util.HichUtil
import hich.rules.ValidationRules
import hich.rules.TransformationRules
import hich.plans.Plan
import hich.plans.Planner

class Engine {
    /* Fast-fail in case HichSpec doesn't have every member of HichSpec.SpecKey
       correctly defined.
    */
    static Map getSpec(SpecKey key) {
        def SPECS = HichSpec.SPECS
        if (!SPECS.containsKey(key)) {
            throw new Exception("FATAL: HichSpec.SPECS does not contain key '${key}'")
        }
        def spec = HichSpec.SPECS[key]
        if (!Map.isInstance(spec)) {
            throw new Exception("FATAL: HichSpec.SPECS[${spec}] is not a Map.")
        }
        return spec
    }

    static Map prepSample(SpecKey key, Map sample) {
        try {
            Map spec = getSpec(key)
            def errors = []

            def PRE_VALIDATORS = PlanKey.PRE_VALIDATORS
            if (spec.containsKey(PRE_VALIDATORS)) {
                errors = validate(spec[PRE_VALIDATORS], sample)
            }
            if (!errors.isEmpty()) {
                throw new Exception("Pre-validation failed on process '${key}': ${errors.join('\n')}\nSample: ${sample}")
            }
            def TRANSFORMERS = PlanKey.TRANSFORMERS
            if (spec.containsKey(TRANSFORMERS)) {
                sample = transform(spec[TRANSFORMERS], sample)
            }

            def POST_VALIDATORS = PlanKey.POST_VALIDATORS
            if (spec.containsKey(POST_VALIDATORS)) {
                errors += validate(spec[POST_VALIDATORS], sample)
            }
            if (!errors.isEmpty()) {
                throw new Exception("Post-validation failed on process '${key}': ${errors.join('\n')}")
            }
        } catch(Exception e) {
            def errors = "Failed to process spec '${key}' on sample ${sample}.\n${e}"
            throw new Exception(errors)
        }
        return sample
    }

    static Map formatOutput(SpecKey key, List outputs) {
        Map spec = getSpec(key)
        def OUTPUTS = PlanKey.OUTPUTS
        
        if (!spec.containsKey(OUTPUTS)) {
            throw new Exception("No OUTPUTS defined in HichSchema for process '${key}'.")
        }
        else {
            def schema = spec[OUTPUTS]
            def count = schema.size()
            def actual_count = outputs.size()
            if (count != actual_count) {
                throw new Exception("In HichSchema, ${key}.OUTPUTS expects ${count} outputs, got ${actual_count}.")
            }
            def result = [:]
            for (int i = 0; i < count; i++) {
                def output_i = outputs[i]
                def schema_i = schema[i]
                if (!schema_i.containsKey("key")) {
                    throw new Exception("In HichSchema, ${key}.OUTPUTS contains schema ${schema_i} which is missing 'key'.")
                }
                def sample_key = schema_i.key
                def sample_value = HichUtil.convertToSchema(output_i, schema_i)
                result += [(sample_key) : sample_value]
            }
            return result
        }
    }

    static List formatInput(SpecKey key, Map sample) {
        def spec = getSpec(key)
        def INPUTS = PlanKey.INPUTS

        if (!spec.containsKey(INPUTS)) {
            throw new Exception("No INPUTS defined in HichSpec for process '${key}': ${spec}.")
        }
        else {
            def inputs = spec[INPUTS]
            def fromSample = inputs.findAll{!it.containsKey("fromProcess") || !it.fromProcess}
            return fromSample.collect {
                input ->
                HichUtil.extractNextflowProcessInput(sample, input)
            }
        }
    }

    static List getInput(SpecKey key, Map sample) {
        sample = prepSample(key, sample)
        return formatInput(key, sample)
    }

    static Plan getPlan(SpecKey key, List input) {
        Map spec = getSpec(key)

        def COMMANDS = PlanKey.COMMANDS
        def INPUTS = PlanKey.INPUTS

        if (!spec.containsKey(COMMANDS)) {
            throw new Error("HichSpec.SPECS.${key}.COMMANDS does not exist.")
        }
        if (!spec.containsKey(INPUTS)) {
            throw new Error("HichSpec.SPECS.${key}.INPUTS does not exist.")
        }

        def planner = spec[COMMANDS].method
        def input_schema = spec[INPUTS]

        def expected_input_arg_types = input_schema.collect{
            if (!it.type) {
                throw new Exception("In SPECS[${key}].INPUTS, required 'type' key not specified for ${it}")
            }
            it.type
        }
        def expected_arg_types = [Map] + expected_input_arg_types

        def error_message
        if (!Planner.metaClass.getMetaMethod(planner, expected_arg_types as Class[])) {
            error_message = "FATAL: Mismatch between HichSpec schema for SPECS[${key}].INPUTS and Planner.${planner} definition.\n"
            
        }

        def commands = spec[COMMANDS].commands
        def builder_input = [commands, *input]
    
        try {
            return Planner."$planner"(*builder_input)
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
                    all_errors += ValidationRules."$rule"(sample, params)
                } catch (Exception e) {
                    def expected_arg_types = [Map, Map] as Class[]
                    def params_class = params.getClass()
                    def actual_arg_types = [Map, params_class]
                    def explanation
                    if (!ValidationRules.metaClass.getMetaMethod(rule, expected_arg_types)) {
                        explanation = "FATAL: HichContract contains misnamed rule '${rule}'.\n"
                    } else if (!ValidationRules.metaClass.getMetaMethod(rule, actual_arg_types)) {
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
                    all_errors += TransformationRules."$rule"(sample, params)
                } catch (Exception e) {
                    def expected_arg_types = [Map, Map] as Class[]
                    def params_class = params.getClass()
                    def actual_arg_types = [Map, params_class]
                    def explanation
                    if (!TransformationRules.metaClass.getMetaMethod(rule, expected_arg_types)) {
                        explanation = "FATAL: HichContract contains misnamed rule '${rule}'.\n"
                    } else if (!TransformationRules.metaClass.getMetaMethod(rule, actual_arg_types)) {
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