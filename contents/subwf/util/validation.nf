def checkMapTypes(Map options, Map validOptions) {
    options.each {
        key, value ->
        // 1. Check if the key is valid
        assert validOptions.containsKey(key), "Invalid key '${key}' found in options. Valid options are: ${validOptions.keySet()}"
        
        // 2. Get the required type from our definition map
        def requiredType = validOptions[key]
        
        // 3. Assert that the value is an instance of the required type
        assert requiredType.isInstance(value), "Option key '${key}' must be of type ${requiredType.simpleName}, but was ${value.getClass().simpleName}"
    }
}