include {formatFromExtension} from './files.nf'

def parsePattern(String str, String parsePattern) {
    /*
        Used to extract sample attributes from filenames, such as condition, biorep, and techrep,
        via a syntax similar to that offered by Python's parse library. In parse, users can
        extract substrings into a map with patterns like: "{condition}_{biorep}_{techrep}.fastq",
        which would take a string like "cond1_br1_tr1.fastq" and return ["condition": "cond1", "biorep": "br1", "techrep": "tr1"].
        This is easier to specify at the command line than a regex but AFAIK has no Groovy equivalent.

        This function implements this parsing functionality, returning the extracted map.
    */

    def patternPlaceholders = []

    // This regex searches the parsePattern string (i.e. "{condition}_{biorep}_{techrep}.fastq")
    // for the placeholders between braces (i.e. condition, biorep, techrep) and adds them
    // to the list of patternPlaceholders to become keys in the output map.

    // It also yields in "pattern" the list of matchers to look for in the input string "str"
    def pattern = parsePattern.replaceAll(/\{([^{}]*)\}/) { match ->
        if (match[1].trim()) {
            patternPlaceholders << match[1]  // Track the placeholder name
            "(?<${match[1]}>.+?)"
        } else {
            ""  // Ignore empty placeholders
        }
    }
    
    // This extracts the patterns from str
    def matcher = str =~ pattern

    // Combine the patternPlaceholders with the corresponding matches from "str" into an output map "result"
    def result = [:]
    if (matcher) {
        patternPlaceholders.eachWithIndex { placeholder, index ->
            result[placeholder] = matcher.group(index + 1)  // Retrieve group by index
        }
    }
    
    return result ?: null
}

def buildFlags(flagsMap) {
    /* Format CLI tool flags passed as hashmap
        
        flagsMap: map of [argName: argVal] pairs
        returns: sing string of combined, single-quoted, space-separated flags.

        Ignores null and false values for flags.
        Converts flags whose value is true to boolean flags.
    */
    flagsMap = flagsMap ?: [:]
    flagsMap = flagsMap.findAll{
        argName, argVal -> 
        def nullVal = argVal == null
        def falseVal = (argVal instanceof Boolean && !argVal)
        def keepVal = (!nullVal) && (!falseVal)
        keepVal
    }
    flagsMap = flagsMap.collect{
        argName, argVal ->
        if (argVal.getClass() == Boolean) {
            argName
        } else {
            "${argName} '${argVal}'"
        }
    }
    flagsMap = flagsMap.join(" ")
    return flagsMap
}

def formatCLIArgs(defaultFlags, updateFlags) {
    def flags = defaultFlags.clone() + (updateFlags ?: [:])
    buildFlags(flags)
}


def formatArg(pattern, object, sep) {
    /*
        Some processes receive either a list of values or a single non-list element
        as parameter values, but need to call a CLI command passing a delimiter-separated
        list of the received values. Other times they get nothing and should
        not pass an argument for that parameter at all. This facilitates this interconversion
        and returns an empty string if the object passed was falsey.

        NOTE this is a potential issue if the goal is to pass a boolean false
        to the CLI command, but I don't think Hich currently does this...

        pattern -- the string pattern to format the results into, like "--numbers {commaSeparatedNums}"
        object -- the element or list of elements to join (where necessary) into a delimiter-separated list
        sep -- the delimiter, like ","
    */
    // Put non-lists into a list so when join is called, it has something to (silently) operate on
    def listed = (object instanceof List || object == null) ? object : [object]
    def joined = listed ? listed.join(sep) : listed

    // Format the string with the result or return 
    return joined ? String.format(pattern, joined) : ""
}

/*
    emptyOnLastStep is used to control behavior when the workflow reaches the step specified by the --lastStep
    param, if specified. It should get passed the main samples channel and returns an empty channel if it's the 
    last step (halting execution) or the original samples channel otherwise.

    If --viewLastStep is specified it will display the contents of the samples channel after the last step finishes.
    If a space-separated string is specified, it will turn that into a list and submap just those sample attributes for viewing.
*/
def emptyOnLastStep(step, samples) {
    def isExplicitLastStep = (params.containsKey("lastStep") && params.get("lastStep") == step)
    def isLastStep = (step == "end") || isExplicitLastStep
    def hasViewLastStep = params.containsKey("viewLastStep") && params.get("viewLastStep")
    if (isLastStep && hasViewLastStep) {
        samples
            | map {
                sample ->
                params.viewLastStep instanceof Boolean ? sample : sample.subMap(params.viewLastStep.split())}
            | view
    }
    return isExplicitLastStep ? channel.empty() : samples
}

def skip(step) {
    /*
        Users may want to skip some steps, such as QC or forming a particular kind of contact matrix,
        or run only certain steps. This uses both params to define a list of steps to be skipped
        (the intersection of skip and runOnly's complement).
    */
    def excluded = params.containsKey("runOnly") && !params.runOnly.split().contains(step)
    def skipped = params.containsKey("skip") && params.skip.split().contains(step)
    return excluded || skipped
}

def sampleFromFastqPairs(files) {
    file1 = files[0]
    file2 = files[1]

    bothFastq = formatFromExtension(file1) == "fastq" && formatFromExtension(file2) == "fastq"
    assert bothFastq, "--fastqPairs grouped ${file1} and ${file2}, but these do not contain a .fastq or .fq extension as expected."
    sample = [fastq1:file1, fastq2: file2, datatype: "fastq"]

    if (params.containsKey("paramsFromPath")) {
        f1Params = parsePattern(file1.toString(), params.paramsFromPath)
        f2Params = parsePattern(file2.toString(), params.paramsFromPath)
        sameParams = f1Params == f2Params
        assert sameParams, "--paramsFromPath yielded different params for ${file1} and ${file2}:\n${f1Params}\n${f2Params}"
        sample += f1Params
    }
    else if (params.containsKey("paramsFromFilename")) {
        f1Params = parsePattern(file1.name.toString(), params.paramsFromFilename)
        f2Params = parsePattern(file2.name.toString(), params.paramsFromFilename)
        sameParams = f1Params == f2Params
        assert sameParams, "--paramsFromFilename yielded different params for ${file1} and ${file2}:\n${f1Params}\n${f2Params}"
        sample += f1Params
    }
    return sample
}