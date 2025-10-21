class HichRuleValidator {
    /**
     * Enforces the parameter contract for a rule or transformer.
     * Throws a fatal exception on failure.
     */
    static void validateParams(String callerName, Map params, Map expectedParams) {
        expectedParams.each { key, expectedType ->
            if (!params.containsKey(key)) {
                throw new Exception("FATAL: HichContract definition for '${callerName}' is missing required parameter '${key}'.")
            }

            def value = params[key]
            if (!(expectedType.isInstance(value))) {
                throw new Exception("FATAL: HichContract parameter '${key}' for '${callerName}' must be type ${expectedType}, but got ${value.class.simpleName}.")
            }
        }
    }
}