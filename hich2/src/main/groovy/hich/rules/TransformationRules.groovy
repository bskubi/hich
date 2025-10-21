package hich.rules
import static hich.util.HichUtil.checkSample

class TransformationRules {
    /* Join values at 'input_keys' by 'separator'
       Unchanged if no 'input_keys' found in 'sample'
    */
    static Map joinPermissive(Map sample, Map params) {
        checkSample(
            params,
            [
                output_key: String,
                input_keys: List,
                separator: String
            ]
        )
        def output_key = params.output_key
        def input_keys = params.input_keys
        def separator = params.separator
        
        if (!input_keys instanceof List) {
            return sample
        }

        def value = (
            input_keys
            .findAll { sample.containsKey(it) }
            .collect { sample[it].toString() }
        )

        if (!value.isEmpty()) {
            sample[output_key] = value.join(separator)
        }

        return sample
    }
}