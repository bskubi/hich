import groovy.transform.Internal
import junit.framework.AssertionFailedError

class CommandChecker {

    /**
     * Splits a command string by the pipe character.
     * @param pipeCmd The full command string containing one or more pipes.
     * @return A list of the separated command parts.
     */
    static List<String> getPipeParts(String pipeCmd) {
        return pipeCmd.split("\\|").collect { it.strip() }
    }

    /**
     * Asserts that a command string starts with an expected name and contains the expected options.
     * @param fullCmdStr The actual command string to check.
     * @param expectedName The expected command name (e.g., "samtools view").
     * @param expectedOpts A list of expected options. Flags are strings (e.g., "-b"),
     * and options with values are maps (e.g., [o: "'test.bam'"]).
     */
    static void checkCmd(String fullCmdStr, String expectedName, List expectedOpts) {
        assert fullCmdStr.startsWith(expectedName), "Command string does not start with '${expectedName}'"
        def actualOpts = getCmdOpts(fullCmdStr, expectedName)
        _checkOpts(actualOpts, expectedOpts)
    }

    /**
     * Extracts options/arguments from a command string by removing the command name.
     * @param fullCmd The command string (e.g., "bwa mem -t 8 file.fq").
     * @param cmdName The command name to remove (e.g., "bwa mem").
     * @return A list of options and arguments.
     */
    static List<String> getCmdOpts(String fullCmd, String cmdName) {
        // Use replaceFirst and trim for robustness
        return fullCmd.replaceFirst(cmdName, "").trim().split(/\s+/).findAll { it }
    }

    /**
     * Internal helper to check for expected options in a list of actual options.
     */
    @Internal // Hides this method from documentation
    private static void _checkOpts(List<String> actualOpts, List expectedItems) {
        expectedItems.each { expected ->
            if (expected instanceof Map) {
                expected.each { key, value ->
                    def actualValue = _getOptVal(actualOpts, key)
                    assert actualValue == value, "For option '${key}', expected value '${value}' but got '${actualValue}'."
                }
            } else { // It's a simple flag
                assert actualOpts.contains(expected), "Expected flag '${expected}' was not found in actual options: ${actualOpts}"
            }
        }
    }

    /**
     * Internal helper to retrieve the value for a given option key (e.g., gets '8' for key '-t').
     */
    @Internal
    private static String _getOptVal(List<String> opts, String name) {
        def index = opts.indexOf(name)
        if (index == -1) {
            throw new AssertionFailedError("Option '${name}' not found in actual options: ${opts}")
        }
        if (index + 1 >= opts.size()) {
            throw new AssertionFailedError("Option '${name}' was found, but it has no value following it in the list.")
        }
        return opts[index + 1]
    }
}