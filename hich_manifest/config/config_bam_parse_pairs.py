from pydantic import FilePath, Field
from .command import Command
from .config_process import ConfigProcess

class ConfigBamParsePairs(ConfigProcess):
    chromsizes: FilePath

    command_samtools_view: Command = Command(base_command="samtools view", flags=["-b"], args=["${bam}"])
    command_pairtools_parse2: Command = Command(base_command="pairtools parse2", opts={"--chroms-path": "${chromsizes}"})
    command_pairtools_sort: Command = Command(base_command="pairtools sort", opts={"--tmpdir": ".", "--nproc-in": "${cpus}", "--nproc-out": "${cpus}", "--output": "${pairs}"})
