import java.nio.file.Path
import java.nio.file.Files


def isExistingFile(it) {
    // Type-agnostic way to check if file exists for any file class having an exists() method.
    return it && it.metaClass.respondsTo(it, 'exists') && it.exists()
}

def formatFromExtension(path) {
    /*
        Look for various known extensions to extract the datatype implicitly from the
        input file so that Hich can ingest intermediate file formats appropriately
        without explicit specification by the user. This is especially helpful in
        permitting the user to use globs at the command line to feed files into Hich.
    */
    def extensions = [".fastq": "fastq",
                  ".fq": "fastq",
                  ".sam": "sambam",
                  ".bam": "sambam",
                  ".pairs": "pairs",
                  ".mcool": "mcool",
                  ".hic": "hic"]
    def pathString = path.toString()
    def foundExtension = extensions.keySet().find {
        ext ->
        pathString.endsWith(ext) || pathString.contains("${ext}.")
    }
    return foundExtension ? extensions[foundExtension] : null
}

def smartFileFinder(pathString) {
    def givenPath = Path.of(pathString)
    def inProjectDir = Path.of("${projectDir}").resolve(pathString)
    return Files.exists(givenPath) ? file(pathString) : file(inProjectDir.toString())
}