package hich.rules
import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.Paths
import static hich.util.HichUtil.checkSample

class ValidationRules {
    static List requireKey(Map map, Map params) {
        List all_errors = checkSample(params, ["key": String])
        
        def key = params.key
        if (!map.containsKey(key)) {
            all_errors += ["Required key '${key}' is missing."]
        }
        return all_errors
    }

    static List notNull(Map map, Map params) {
        List all_errors = checkSample(params, ["key": String])

        def key = params.key
        if (map.containsKey(key) && map[key] == null) {
            all_errors += ["Value at key '${key}' must not be null."]
        }
        return all_errors
    }

    static List classChoice(Map map, Map params) {
        List all_errors = checkSample(
            params, 
            [
                "key": String,
                "classes": List<Class>
            ]
        )

        def key = params.key
        def classes = params.classes

        if (!map.containsKey(key)) {
            return []
        }

        def value = map[key]

        if (value == null) {
            return []
        }

        def type_match = classes.any{
            cls ->
            cls.isInstance(value)
        }
        
        if (!type_match) {
            def classSimpleNames = classes.collect{it.simpleName}
            String error_message = (
                "Value at key '${key}' is disallowed class " +
                "${value.class.simpleName}. Must be one of the following " +
                "classes: ${classSimpleNames}."
            )
            all_errors += [error_message]
        }

        return all_errors
    }

    static List hasNonWhitespaceChar(Map map, Map params) {
        List all_errors = checkSample(params, ["key": String])

        def key = params.key

        if (!map.containsKey(key)) {
            return []
        }

        def value = map[key]

        if (!String.isInstance(value)) {
            return []
        }

        if (value.strip() == "") {
            all_errors += ["Value at key '${key}' does not contain non-whitespace chars."]
        }

        return all_errors
    }

    static List isIn(Map map, Map params) {
        List all_errors = checkSample( 
            params, 
            [
                "key": String,
                "choices": List
            ]
        )

        def key = params.key
        if (!map.containsKey(key)) {
            return []
        }

        def value = map[key]
        def choices = params.choices
        if (!choices.contains(map[key])) {
            all_errors += ["Value at key '${key}' is not one of: ${choices}"]
        }
        
        return all_errors
    }

    static List pathCompatibleClass(Map map, Map params) {
        List all_errors = checkSample(params, ["key": String])
        
        def key = params.key
        if (!map.containsKey(key)) {
            return []
        }

        def value = map[key]

        if (!Path.isInstance(value) && !String.isInstance(value)) {
            String error_message = (
                "Value at key '${key}' must be a Path or String, but is an " +
                "${value.class.simpleName}"
            )
            all_errors += [error_message]
        }

        return all_errors
    }

    static List pathExists(Map map, Map params) {
        List all_errors = checkSample(params, ["key": String])

        def key = params.key
        if (!map.containsKey(key)) {
            return []
        }

        def value = map[key]
        Path path_to_check

        if (Path.isInstance(value)) {
            path_to_check = value
        } else if (String.isInstance(value)) {
            path_to_check = Paths.get(value)
        } else {
            return []
        }

        if (!Files.exists(path_to_check)) {
            def absolute_path = path_to_check.toAbsolutePath()
            String error_message = (
                "Value at key '${key}' points to a path that does not " +
                "exist: ${absolute_path}"
            )
            
            all_errors += [error_message]
        }

        return all_errors
    }

    static List notNullCountRange(Map map, Map params) {
        List all_errors = checkSample(params, ["keys": List<String>, "allowed": List<Integer>])
        def keys = params.keys
        def allowed = params.allowed
        def count = keys.findAll{map.getAt(it) != null}.size()

        if (!allowed.contains(count)) {
            String error_message = (
                "Map contained ${count} keys from set ${keys}. " +
                "Allowed counts: ${allowed}"
            )
            all_errors += [error_message]
        }

        return all_errors
    }

    static List requireExistingPath(Map map, Map params) {
        def all_errors = []
        all_errors += requireKey(map, params)
        all_errors += notNull(map, params)
        all_errors += pathCompatibleClass(map, params)
        all_errors += pathExists(map, params)
        return all_errors
    }

    static List requireIsIn(Map map, Map params) {
        def all_errors = []
        all_errors += requireKey(map, params)
        all_errors += isIn(map, params)
        return all_errors
    }
}