import groovy.text.SimpleTemplateEngine

class HichUtil {
    static Map formatIfPresent(Map pattern_map, String pattern_key, Map bind) {
        try {
            if (pattern_map.containsKey(pattern_key) && pattern_map[pattern_key] != null) {
                def engine = new SimpleTemplateEngine()
                def pattern = pattern_map[pattern_key]
                def formatted = engine.createTemplate(pattern).make(bind).toString()
                pattern_map += [(pattern_key): formatted]
            }
        } catch(Exception e) {
            throw new Exception("FATAL: Could not format with pattern_key: '${pattern_key}', pattern: '${pattern_map[pattern_key]}', bind: ${bind}\n${e}")
        }

        return pattern_map
    }



    static String format(String pattern, Map bind) {
        def engine = new SimpleTemplateEngine()
        return engine.createTemplate(pattern).make(bind).toString()
    }
}