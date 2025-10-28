from typing import Union, Optional
from typing_extensions import Self
from enum import Enum
from pydantic import BaseModel, Field, FilePath, DirectoryPath, model_validator, field_validator, AliasChoices

from .process_config import ProcessConfig
from .command import Command

class BWA(str, Enum):
    MEM = "bwa mem"
    MEM2 = "bwa-mem2 mem"

class BWAMETH(str, Enum):
    MEM = "python -m bwameth"
    MEM2 = "python -m bwameth-mem2"

class FastqAlign(ProcessConfig):
    fastq1: Optional[FilePath] = Field(None)
    fastq2: Optional[FilePath] = Field(None)
    aligner: Optional[Union[BWA,BWAMETH]] = Field(None)
    aligner_index_dir: Optional[DirectoryPath] = Field(None)
    aligner_index_prefix: str = Field(None)

    command_align: Optional[Command] = Field(None)
    command_samtools_view: Optional[Command] = Field(None)

    @model_validator(mode = 'before')
    @classmethod
    def set_defaults(cls, data: dict) -> dict:
        fastq1 = data.get('fastq1')
        fastq2 = data.get('fastq2')
        aligner = data.get('aligner')
        aligner_index_dir = data.get('aligner_index_dir')
        aligner_index_prefix = data.get('aligner_index_prefix')        
        
        if all([fastq1, aligner, aligner_index_dir, aligner_index_prefix]):
            if not data.get('command_align'):
                # Build command
                command_align = {'base_command': aligner, 'opts': {}, 'flags': [], 'args': []}

                # Deal with fastq data files
                if fastq2:
                    command_align['flags'] += ['-p']
                    command_align['args'] = ['${fastq1}', '${fastq2}']
                else:
                    command_align['args'] = ['${fastq1}']
                
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
        elif not (fastq1 or fastq2 or aligner or aligner_index_dir or aligner_index_prefix):
            data['skip'] = True
        else:
            raise ValueError(f"Specified some but not all of fastq1, fastq2, aligner, aligner_index_dir, aligner_index_prefix.")
        return data