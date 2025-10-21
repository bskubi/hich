package hich.util
import java.nio.file.Path
import java.nio.file.Paths
/* Returns list of errors for missing keys, mismatched types.
*/
class HichUtil {
    static List checkSample(Map sample, Map requiredKeysTypes) {
        List all_errors = []

        requiredKeysTypes.each { key, expectedType ->
            // Key must be in sample
            if (!sample.containsKey(key)) {
                all_errors.add(
                    "FATAL: Spec for '${name}' missing required key '${key}'."
                )
                return
            }

            // Value of key must be type
            def value = sample[key]
            if (!(expectedType.isInstance(value))) {
                all_errors.add(
                    "FATAL: Spec key '${key}' " +
                    "must be type ${expectedType}, but got " +
                    "${value.class.simpleName}."
                )
            }
        }

        return all_errors
    }

    /* Convert object to Path type
    */
    static Path toPath(Object path) {
        return Paths.get(path.toString())
    }

    /* Check if Path, String, etc. is existing path
    */
    static Boolean pathExists(Object path) {
        return toPath(path).exists()
    }

    /* 
    */
    static def extractNextflowProcessInput(Map sample, Map input) {
        List all_errors = []
        if (!input.containsKey("key")) {
            all_errors.add(
                "Input does not contain required key 'key'"
            )
        }

        String key = input.key
        Boolean sampleContainsKey = sample.containsKey(key)
        if (!sampleContainsKey && !input.containsKey("whenMissing")) {
            throw new Exception("Sample did not contain required key '${key}'.")
        }

        def value = sampleContainsKey ? sample[key] : input.whenMissing

        if (!input.containsKey("type")) {
            throw new Exception("Input ${input} did not contain required key 'type'.")
        }
        if (input.type == Path) {
            value = sampleContainsKey ? toPath(value) : value
        }

        return value
    }

    static def convertToSchema(item, Map schema) {
        def result = item
        if (schema.containsKey("type")) {
            def itemType = schema.type
            if (!Class.isInstance(itemType)) {
                throw new Exception("Schema ${schema} contains 'type' ${itemType} that doesn't store a Class")
            }
            result = item.asType(itemType)
        }
        return result
    }


}
