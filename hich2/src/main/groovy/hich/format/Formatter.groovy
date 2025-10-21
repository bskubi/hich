package hich.format
import groovy.text.SimpleTemplateEngine
/* Format raw pattern String containing GString, i.e. '${aligner}'
*/
class Formatter {
    String format(String pattern, Map bind) {
        def engine = new SimpleTemplateEngine()
        return engine.createTemplate(pattern).make(bind).toString()
    }

    /* Return 'patterns' after formatting 'key' with 'bind'.
    */
    Map updateFormat(Map patterns, String key, Map bind) {
        try {
            Boolean hasPattern = (
                patterns.containsKey(key) &&
                patterns[key] != null
            )
            if (hasPattern) {
                String pattern = patterns[key]
                String formatted = format(pattern, bind)
                patterns += [(key): formatted]
            }
        } catch(Exception e) {
            String error = (
                "FATAL: Could not update format for " +
                "key: '${key}', " +
                "pattern: '${patterns[key]}', " +
                "bind: ${bind}\n" +
                "${e}"
            )
            throw new Exception(error)
        }

        return patterns
    }
}




