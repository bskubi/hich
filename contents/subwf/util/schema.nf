import groovy.yaml.YamlSlurper
import groovy.json.JsonOutput

def validation_error(sample, name, message, schema) {
    def error_message = [
        "Attribute": name,
        "Error": name + message,
        "Schema": schema,
        "Sample:": sample
    ].findAll{k, v -> v}

    def jsonError = JsonOutput.toJson(error_message)
    jsonError = JsonOutput.prettyPrint(jsonError)
    error(jsonError)
}

def validateType(value, okTypes) {
    def types = [
        "String": String,
        "List": List,
        "Map": Map,
        "Integer": Integer,
        "Float": Float,
        "Path": Path
    ]
    if (okTypes instanceof String) {
        okTypes = [okTypes]
    }
    if (value == null) {
        return okTypes.contains("Null")
    }
    
    // Use collect to transform and chain any directly
    return okTypes
        .collect { types[it] }
        .findAll()
        .any {
            expectedClass ->
            if (expectedClass == Path) {
                if (value instanceof Path) {
                    return value.exists()
                } else {
                    return file(value).exists()
                }
            } 
            return expectedClass.isInstance(value)
        }
}

def validate(schemas, sample) {
    def validated = [:]
    schemas.each {
        name, schema ->
        def value = sample[name]
        
        if (!schema.missingOK && !sample.containsKey(name) && !schema.containsKey("default")) {
            validation_error(sample, name, " is missing (required, no default)", schema)
        }
        if (!sample.containsKey(name) && schema.containsKey("default")) {
            value = schema.default
        }
        if (!sample.containsKey(name) && schema.missingOK) {
            return
        }
        if (schema.type && !validateType(value, schema.type)) {
            if (schema.type.contains("Path")) {
                def absPath = value == null ? null : file(value).toAbsolutePath()
                def message = " is type ${value.getClass()}, but should be one of: ${schema.type}. As a path, resolves to ${absPath}, which does not exist."
            } else {
                def message = " is type ${value.getClass()}, but should be one of: ${schema.type}"
            }
            validation_error(sample, name, message, schema)
        }
        if (schema.ignoreNull && sample.containsKey(name) && sample[name] == null) {
            return
        }
        if (schema.min != null && value < schema.min) {
            validation_error(sample, name, " = ${value} < ${schema.min} (min)", schema)
        }
        if (schema.max != null && value > schema.max) {
            validation_error(sample, name, " = ${value} > ${schema.max} (max)", schema)
        }
        if (schema.choices != null && !schema.choices.contains(value)) {
            validation_error(sample, name, " = ${value} not in ${schema.choices}", schema)
        }
        if (sample.containsKey(name)) {
            if (schema.require_attribs) {
                schema.require_attribs.each {
                    required ->
                    if (!sample.containsKey(required)) {
                        validation_error(sample, name, " requires ${required} missing from sample", schema)
                    }
                }
            }
            if (schema.forbid_attribs) {
                schema.forbid_attribs.each {
                    forbidden ->
                    if (sample.containsKey(forbidden)) {
                        validation_error(sample, name, " forbids ${forbidden} present in sample", schema)
                    }
                }
            }
        }

        validated += [(name): value]
    }
    return validated
}

def extract_inputs(process, config, sample) {
    def input = config[process].input
    def schemas = input.collectEntries{
        name, attribs ->
        // Get schema from source group if specified and
        // override source schema with local definitions
        def source_schema = attribs.source ? config[attribs.source][name] : [:]
        def schema = (source_schema ?: [:]) + attribs.findAll{k, v -> k != "source"}
        [name, schema]
    }
    return validate(schemas, sample)
}

def load_sample_attrib_schema() {
    def sample_attrib_schema_file = File(params.get("sampleAttribSchema", "sample-attrib-schema.yaml"))
    def sample_attrib_schema = new YamlSlurper().parseText(sample_attrib_schema_file.text)
    return sample_attrib_schema
}