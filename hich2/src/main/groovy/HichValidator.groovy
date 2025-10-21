import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.Paths

class HichValidator {
    static List requireKey(Map map, Map params) {
        HichRuleValidator.validateParams("requireKey", params, ["key": String])
        
        def key = params.key
        if (!map.containsKey(key)) {
            return ["Required key '${key}' is missing."]
        }
        return []
    }

    static List notNull(Map map, Map params) {
        HichRuleValidator.validateParams("notNull", params, ["key": String])

        def key = params.key
        if (map.containsKey(key) && map[key] == null) {
            return ["Value at key '${key}' must not be null."]
        }
        return []
    }

    static List classChoice(Map map, Map params) {
        HichRuleValidator.validateParams(
            "classChoice", 
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
            return ["Value at key '${key}' is disallowed class ${value.class.simpleName}. Must be one of the following classes: ${classSimpleNames}."]
        }

        return []
    }

    static List hasNonWhitespaceChar(Map map, Map params) {
        HichRuleValidator.validateParams("hasNonWhitespaceChar", params, ["key": String])

        def key = params.key

        if (!map.containsKey(key)) {
            return []
        }

        def value = map[key]

        if (!String.isInstance(value)) {
            return []
        }

        if (value.strip() == "") {
            return ["Value at key '${key}' does not contain non-whitespace chars."]
        }

        return []
    }

    static List isIn(Map map, Map params) {
        HichRuleValidator.validateParams(
            "isIn", 
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
            return ["Value at key '${key}' is not one of: ${choices}"]
        }
        
        return []
    }

    static List pathCompatibleClass(Map map, Map params) {
        HichRuleValidator.validateParams("pathCompatibleClass", params, ["key": String])
        
        def key = params.key
        if (!map.containsKey(key)) {
            return []
        }

        def value = map[key]

        if (!Path.isInstance(value) && !String.isInstance(value)) {
            return ["Value at key '${key}' must be a Path or String, but is an ${value.class.simpleName}"]
        }

        return []
    }

    static List pathExists(Map map, Map params) {
        HichRuleValidator.validateParams("pathExists", params, ["key": String])

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
            return ["Value at key '${key}' points to a path that does not exist: ${absolute_path}"]
        }

        return []
    }

    static List notNullCountRange(Map map, Map params) {
        HichRuleValidator.validateParams("pathExists", params, ["keys": List<String>, "allowed": List<Integer>])
        def all_errors = []
        def keys = params.keys
        def allowed = params.allowed
        def count = keys.findAll{map.getAt(it) != null}.size()

        if (!allowed.contains(count)) {
            return ["Map contained ${count} keys from set ${keys}. Allowed counts: ${allowed}"]
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