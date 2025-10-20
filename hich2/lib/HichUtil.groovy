import groovy.text.SimpleTemplateEngine

class HichUtil {
    static Map formatIfPresent(Map pattern_map, String pattern_key, Map bind) {
        if (pattern_map.containsKey(pattern_key) && pattern_map[pattern_key] != null) {
            def engine = new SimpleTemplateEngine()
            def pattern = pattern_map[pattern_key]
            print(pattern)
            def formatted = engine.createTemplate(pattern).make(bind).toString()
            print(formatted)
            pattern_map += [(pattern_key): formatted]
        }
        return pattern_map
    }

    static String format(String pattern, Map bind) {
        def engine = new SimpleTemplateEngine()
        return engine.createTemplate(pattern).make(bind).toString()
    }
}