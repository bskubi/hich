``util/cli.nf``
--------------------------

.. code-block:: c

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