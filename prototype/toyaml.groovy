#!/usr/bin/env groovy
@Grab(group='org.yaml', module='snakeyaml', version='1.29')
import org.yaml.snakeyaml.Yaml
import groovy.util.ConfigSlurper

// Function to load a Nextflow/Groovy-style config file as a Groovy object
def loadConfig(String filePath) {
    // Parse the config file using ConfigSlurper
    def config = new ConfigSlurper().parse(new File(filePath).toURI().toURL())
    return config
}

// Function to convert a Groovy object to YAML and print to stdout
def convertToYaml(def groovyObject) {
    Yaml yaml = new Yaml()
    String yamlOutput = yaml.dump(groovyObject)  // Convert Groovy object to YAML format
    println(yamlOutput)
}

// Main logic to accept a file from the CLI and process it
if (args.length != 1) {
    println "Usage: groovy toyaml.groovy <input_file>"
    System.exit(1)
}

String filePath = args[0]
def configObject = loadConfig(filePath)  // Load the config file as a Groovy object
convertToYaml(configObject)  // Convert and print the object as YAML
