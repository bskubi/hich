import java.nio.file.Path
import java.lang.reflect.Modifier

class HichContract {
    enum Schema {
        PRE_VALIDATORS, 
        TRANSFORMERS, 
        POST_VALIDATORS, 
        INPUTS, 
        COMMANDS,
        OUTPUTS
    }

    static {
        // Get all enum values, convert them to strings.
        def valid_keys = Schema.values().collect { it.toString() }.toSet()

        if (!PROCESS_INTERFACES) {
            throw new Exception("FATAL: Contract.PROCESS_INTERFACES does not exist.")
        }
        if (!Map.isInstance(PROCESS_INTERFACES)) {
            throw new Exception("FATAL: Contract.PROCESS_INTERFACES is not a Map.")
        }

        PROCESS_INTERFACES.each { processName, interfaceDef ->
            // The map keys are enums, so we convert them to strings for the check.
            def present_keys = interfaceDef.keySet().collect { it.toString() }.toSet()
            def unknown_keys = present_keys - valid_keys
            if (!unknown_keys.isEmpty()) {
                throw new Exception("FATAL: HichContract for process '${processName}' contains unknown keys: ${unknown_keys}")
            }
        }
    }

    private static final Map ATTRIBUTES = [
        id:                 [type: String,  doc: "Sample ID"],
        condition:          [type: String,  doc: "Experimental condition ID"],
        biorep:             [type: String,  doc: "Biological replicate ID"],
        assembly:           [type: String,  doc: "Name of reference genome."],
    ]

    private static final Map PROCESS_INTERFACES = [
        ADD_SAMPLE : [
            (Schema.TRANSFORMERS): [
                joinPermissive: [
                    [
                        output_key: "id", 
                        input_keys: ["condition", "biorep", "techrep"], 
                        separator: "_"
                    ]
                ]
            ],
            (Schema.POST_VALIDATORS): [
                hasNonWhitespaceChar: [
                    [key: "id"],
                    [key: "assembly"]
                ]
            ]
        ],
        ALIGN : [
            (Schema.PRE_VALIDATORS): [
                requireExistingPath: [ [key: "aligner_index_dir"] ],
                requireIsIn: [ [key: "aligner", choices: ["bwa mem", "bwa-mem2", "bwameth", "bwameth-mem2"]] ],
                requireKey: [[key: "aligner_index_prefix"]],
                notNull: [ [key: "aligner_opts"], [key: "aligner_index_prefix"] ],
                classChoice: [ [key: "aligner_opts", classes: [Map]] ],
                notNullCountRange: [ 
                    [ keys: ["aligner_opts.bwa", "aligner_opts.bwameth"], allowed: [0, 1] ],
                    [ keys: ["fastq", "fastq1"], allowed: [1] ],
                    [ keys: ["fastq1", "fastq2"], allowed: [0, 2] ]
                ]
            ],
            (Schema.INPUTS): [
                [key: "id", type: String],
                [key: "fastq", type: Path, whenMissing: []],
                [key: "fastq1", type: Path, whenMissing: []],
                [key: "fastq2", type: Path, whenMissing: []],
                [key: "aligner", type: String],
                [key: "aligner_index_dir", type: Path],
                [key: "aligner_index_prefix", type: String],
                [key: "aligner_opts", type: Map, whenMissing: null],
                [key: "cpus", type: Integer, fromProcess: true]
            ],
            (Schema.COMMANDS): [
                method: "alignContext",
                commands: [
                    bwa_mem_command: [
                        base_command: '${aligner}',
                        default_options: [
                            "-S": true,
                            "-P": true,
                            "-5": true,
                            "-M": true,
                            "-p": 'fastq != null',
                            "-t": '${cpus}'
                        ],
                        arguments: [
                            "aligner_index": '${aligner_index_dir}/${aligner_index_prefix}',
                            "fastq": 'fastq ?: null',
                            "fastq1": 'fastq1 ?: null',
                            "fastq2": 'fastq2 ?: null'
                        ]
                    ],
                    bwameth_command: [
                        base_command: '${aligner}',
                        default_options: [
                            "--do-not-penalize-chimeras": true,
                            "--reference": '${aligner_index_dir}/${aligner_index_prefix}',
                            "-p": 'fastq != null',
                            "-t": '${cpus}'
                        ],
                        arguments: [
                            "fastq": 'fastq ?: null',
                            "fastq1": 'fastq1 ?: null',
                            "fastq2": 'fastq2 ?: null'
                        ]
                    ],
                    samtools_view_command: [
                        base_command: "samtools view",
                        required_options: [
                            "-b": true,
                            "-o": '${id}.bam'
                        ]
                    ]
                ]
            ],
            (Schema.OUTPUTS): [
                [key: "id", type: String],
                [key: "bam", type: String]
            ]
        ]
    ]
}