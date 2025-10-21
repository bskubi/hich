package hich.specs.plans
import hich.specs.plans.keys.PlanKey
import java.nio.file.Path

class Align {
    private static final Map spec = [
        (PlanKey.PRE_VALIDATORS): [
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
        (PlanKey.INPUTS): [
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
        (PlanKey.COMMANDS): [
            method: "planAlign",
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
        (PlanKey.OUTPUTS): [
            [key: "id", type: String],
            [key: "bam", type: String]
        ]
    ]
}