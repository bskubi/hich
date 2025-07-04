def validateMap(Map options, Map validOptions, List<String> checks) {
    if (checks.contains("requireKeys")) {
        // 1. Check that required keys are present
        validOptions.keySet().each {
            reqKey ->
            assert options.containsKey(reqKey), "Missing required key '${reqKey}' from options ${options}."
        }
    }

    options.each {
        key, value ->
        
        if (checks.contains("limitKeys")) {
            // 2. Check if the key is valid
            assert validOptions.containsKey(key), "Invalid key '${key}' found in options. Valid options are: ${validOptions.keySet()}"
        }

        if (checks.contains("limitTypes")) {
            // 3. Assert that the value is an instance of the required type
            def requiredType = validOptions[key]

            if (requiredType != null) {

                if (!(requiredType instanceof List)) {
                    requiredType = [requiredType]
                }

                requiredType.each {
                    rqType ->
                    
                    assert rqType.isInstance(value), "Option key '${key}' must be one of the following types: ${requiredType.simpleName} but was an ${value.getClass().simpleName} ${value}"
                }
            }
        }
    }
}