from typing import Union, Optional
from enum import Enum
from pydantic import model_validator

from .config_process import ConfigProcess
from .command import Command

class BWA(str, Enum):
    MEM = "bwa mem"
    MEM2 = "bwa-mem2 mem"

class BWAMETH(str, Enum):
    MEM = "python -m bwameth"
    MEM2 = "python -m bwameth-mem2"

class FASTQ_TYPE(str, Enum):
    PAIRED_END = "PAIRED_END"
    SINGLE_END = "SINGLE_END"
    INTERLEAVED = "INTERLEAVED"

class ConfigFastqAlign(ConfigProcess):
    fastq_type: FASTQ_TYPE
    aligner: Optional[Union[BWA,BWAMETH]] = None
    aligner_index_prefix: str = None

    command_align: Optional[Command] = None
    command_samtools_view: Optional[Command] = None

    @model_validator(mode = 'before')
    @classmethod
    def set_defaults(cls, data: dict) -> dict:
        fastq = data.get('fastq')
        aligner = data.get('aligner')
        aligner_index_prefix = data.get('aligner_index_prefix')        
        
        if all([fastq, aligner, aligner_index_prefix]):
            if not data.get('command_align'):
                # Build command
                command_align = {'base_command': aligner, 'opts': {}, 'flags': [], 'args': []}

                # Deal with fastq data files
                if fastq == FASTQ.PAIRED_END:
                    command_align['args'] = ['${fastq1}', '${fastq2}']
                else:
                    command_align['args'] = ['${fastq1}']
                if fastq == FASTQ.INTERLEAVED:
                    command_align['flags'] += ['-p']
                
                # Deal with aligner-specific flags
                aligner_index = '${aligner_index_dir}/' + aligner_index_prefix
                if isinstance(aligner, BWA):
                    command_align['flags'] += ["-S", "-P", "-5", "-M"]
                    command_align['args'] = [aligner_index] + command_align['args']
                elif isinstance(aligner, BWAMETH):
                    command_align['opts'] = {"--reference": aligner_index}
                    command_align['flags'] += ["--do-not-penalize-chimeras"]
                
                # Build command
                data['command_align'] = Command(**command_align)

            if not data.get('command_samtools_view'):
                data['command_samtools_view'] = Command(base_command="samtools view", opts={"-o": "${bam}.bam"}, flags=["-b"])
        return data
