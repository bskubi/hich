package hich.errors
import hich.format.Formatter

class Errors {
    String error_pattern = '${error}'
    List errors = []
    Boolean numbered = true

    def AddError(String error) {
        String error_message = Formatter.format(error_pattern, [error: error])
        if (numbered) {
            String N = (errors.size() + 1).toString()
            error_message = "${N}. ${error_message}"
        }
        errors.add(error_message)
    }

    def ThrowIfErrors() {
        if (errors) {
            String error_message = errors.join("\n")
            throw new Exception(errors)
        }
    }
}