import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.Paths

class HichWorkflowAdapter {
    
    static def toPath(String value) {
        def Path path

        if (Path.isInstance(value)) {
            path = value
        } else if (String.isInstance(value)) {
            path = Paths.get(value)
        } else {
            throw new Exception("Could not convert '${value}' to file.")
        }

        if (!Files.exists(path)) {
            def absolute_path = path.toAbsolutePath()
            return ["Value '${value}' points to a path that does not exist: ${absolute_path}"]
        }

        return path
    }

    static def extractNextflowProcessInput(Map sample, Map input) {
        if (!input.containsKey("key")) {
            throw new Exception("Input does not contain required key 'key'")
        }

        def key = input.key
        def sampleContainsKey = sample.containsKey(key)

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

    static def toPath(Object taskPath) {
        return taskPath ? Paths.get(taskPath.toString()) : null
    }
}