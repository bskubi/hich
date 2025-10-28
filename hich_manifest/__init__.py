# from typing import Union, Optional
# from typing_extensions import Self
# from enum import Enum
# from pydantic import BaseModel, Field, FilePath, DirectoryPath, model_validator, field_validator, AliasChoices

# class Command(BaseModel):
#     base_command: Union[Enum, str]
#     opts: dict = Field({})
#     flags: list[str] = Field([])
#     args: list[str] = Field([])

# class Entrypoint(str, Enum):
#     FASTQ_ALIGN = "FASTQ_ALIGN"
#     MERGE_BEFORE_DEDUP = "PAIRS_MERGE_BEFORE_DEDUP"
#     MERGE_AFTER_DEDUP = "PAIRS_MERGE_AFTER_DEDUP"

# class BWA(str, Enum):
#     MEM = "bwa mem"
#     MEM2 = "bwa-mem2 mem"

# class BWAMETH(str, Enum):
#     MEM = "python -m bwameth"
#     MEM2 = "python -m bwameth-mem2"

# class ProcessConfig(BaseModel):
#     skip: bool = False

# class FastqAlign(ProcessConfig):
#     fastq1: Optional[FilePath] = Field(None)
#     fastq2: Optional[FilePath] = Field(None)
#     aligner_index_dir: Optional[DirectoryPath] = Field(None)

#     command_align: Command = None
#     command_samtools_view: Command = Field(Command(base_command="samtools view", opts={"-o": "${bam}.bam"}, flags=["-b"]))

#     @model_validator(mode = 'before')
#     @classmethod
#     def set_command_align(cls, data: dict) -> dict:
#         fastq1 = data.get('fastq1')
#         fastq2 = data.get('fastq2')
#         aligner = data.get('aligner')
#         aligner_index_dir = data.get('aligner_index_dir')
#         aligner_index_prefix = data.get('aligner_index_prefix')
#         if aligner is None:
#             data['command_align'] = 

#         try:
#             aligner = BWA(aligner)
#         except:
#             aligner = BWAMETH(aligner)
        

#         if data.get('command_align') is None and all([fastq1, aligner, aligner_index_dir, aligner_index_prefix]):
#             aligner_index = '${aligner_index_dir}/' + aligner_index_prefix

#             command_align = {"base_command": aligner}
#             if fastq2:
#                 command_align['flags'] = []
#                 command_align['args'] = ['${fastq1}', '${fastq2}']
#             else:
#                 command_align['flags'] = ['-p']
#                 command_align['args'] = ['${fastq1}']

#             if isinstance(aligner, BWA):
#                 command_align['flags'] += ["-S", "-P", "-5", "-M"]
#                 command_align['args'] = [aligner_index] + command_align['args']
#             else:
#                 command_align['opts'] = {"--reference": aligner_index}
#                 command_align['flags'] += ["--do-not-penalize-chimeras"]
            
#             command_align = Command(**command_align)
#         data['command_align'] = command_align
#         return data

# class BamParsePairs(ProcessConfig):
#     chromsizes: FilePath

#     command_samtools_view: Command = Field(Command(base_command="samtools view", flags=["-b"]))
#     command_pairtools_parse2: Command = Field(Command(base_command="pairtools parse2", opts={"--chroms-path": "${chromsizes}"}))
#     command_pairtools_sort: Command = Field(Command(base_command="pairtools sort", opts={"--tmpdir": ".", "--nrpoc-in": "${cpus}", "--nproc-out": "${cpus}", "--output": "${pairs}"}))

# class PairsLabel(ProcessConfig):
#     fragment_index: Optional[FilePath] = Field(None)

#     command_hich_pairs_map_ends: Optional[Command] = Field(None)

#     @model_validator(mode='before')
#     @classmethod
#     def set_defaults(cls, data: dict) -> dict:
#         if data.get('fragment_index') and not data.get('command_hich_pairs_map_ends'):
#             data['command_hich_pairs_map_ends'] = Command(
#                 base_command="hich pairs map-ends",
#                 opts={
#                     "--idx1": "rfrag1",
#                     "--start1": "rfrag_start1",
#                     "--end1": "rfrag_end1",
#                     "--idx2": "rfrag2",
#                     "--start2": "rfrag_start2",
#                     "--end2": "rfrag_end2"
#                 },
#                 args=["${fragment_index}", "${pairs}", "${pairs_labeled}"]
#             )
#         else:
#             data['skip'] = True

    
# class PairsSelect(ProcessConfig):
#     command_pairtools_select: Command = Field(
#         Command(
#             base_command = "pairtools select",
#             opts = {
#                 "--output": "${pairs_selected}",
#                 "--nproc-in": "${cpus}",
#                 "--nproc-out": "${cpus}"
#             },
#             args = ["true", "${pairs}"]
#         )
#     )

# class PairsDedup(ProcessConfig):
#     command_pairtools_dedup: Command = Field(
#         Command(
#             base_command = "pairtools dedup",
#             opts = {
#                 "--output": "${pairs_deduped}",
#                 "--nproc-in": "${cpus}",
#                 "--nproc-out": "${cpus}"
#             },
#             args = ["${pairs}"]
#         )
#     )

# class PairsCoolBin(ProcessConfig):
#     chromsizes: FilePath

#     command_cooler_cload_pairs: Command = Field(
#         Command(
#             base_command = "cooler cload pairs",
#             opts = {
#                 "-c1": "2",
#                 "-p1": "3",
#                 "-c2": "4",
#                 "-p2": "5",
#                 "--temp-dir": "."
#             },
#             args = ["${chromsizes}:1000", "${pairs}", "${cool}"]
#         )
#     )

# class CoolCoarsenAddnorm(ProcessConfig):
#     command_cooler_zoomify: Command = Field(
#         Command(
#             base_command = "cooler zoomify",
#             opts = {
#                 "--nproc": "${cpus}",
#                 "--balance-args": "--nproc ${cpus} --max-iters 2000",
#                 "--resolutions": "1000,2000,5000,10000,20000,50000,100000",
#                 "--out": "${mcool}"
#             },
#             flags = ["--balance"],
#             args = ["${mcool}"]
#         )
#     )

# class PairsHiCBinCoarsenAddnorm(ProcessConfig):
#     command_juicer_tools_pre: Command = Field(
#         Command(
#             base_command = "juicer_tools pre",
#             opts = {
#                 "-r": "1000,2000,5000,10000,20000,50000,100000",
#                 "--threads": "${cpus}",
#                 "-t": "."
#             },
#             flags = ["-Xms${xms}g", "-Xmx${xmx}g"],
#             args = ["${pairs}", "${hic}", "${chromsizes}"]
#         )
#     )

# class PairsMerge(ProcessConfig):
#     merge_ids: list[str] = Field(min_length=1)
    
#     command_pairtools_merge: Command = Field(
#         Command(
#             base_command = "pairtools merge",
#             opts = {
#                 "--output": "${pairs_merged}",
#                 "--nproc-in": "${cpus}",
#                 "--nproc-out": "${cpus}",
#                 "--tmpdir": "."
#             },
#             args = ["${pairs}"]
#         )
#     )

# class Fastq(BaseModel):
#     fastq1: FilePath

# class FastqPairedEnd(Fastq):
#     fastq2: FilePath

# class FastqSingleEnd(Fastq): ...

# class FastqInterleaved(Fastq): ...
    
# class Bam(BaseModel):
#     bam: FilePath

# class Pairs(BaseModel):
#     pairs: FilePath

# class PairsList(BaseModel):
#     pairs: list[Pairs]

# class Cool(BaseModel):
#     cool: FilePath

# class Mcool(BaseModel):
#     mcool: FilePath

# class Hic(BaseModel):
#     hic: FilePath

# def set_fastq_align(data: dict) -> dict:
#     input_data = data.get('input_data')
#     if isinstance(input_data, FastqPairedEnd):
#         data['fastq2'] = input_data.fastq2
#     if isinstance(input_data, Fastq):
#         data['fastq1'] = input_data.fastq1
#         if not data.get('fastq_align'):
#             data['fastq_align'] = FastqAlign(**data)
#     elif data.get('fastq_align') is None:
#         data['fastq_align'] = FastqAlign(skip=True)
#     return data

# class Record(BaseModel):
#     input_data: Union[FastqPairedEnd, FastqSingleEnd, FastqInterleaved, Bam, Pairs, Cool, Mcool, Hic]
#     entrypoint: Entrypoint = Field(None)
#     fastq_align: Optional[FastqAlign] = Field(None)

#     @model_validator(mode = 'before')
#     @classmethod
#     def set_defaults(cls, data: dict) -> dict:
#         return set_fastq_align(data)



# fastq1 = "tests/assets/fastq/1k/1k_ERR1413593_1.fastq.gz"
# fastq2 = "tests/assets/fastq/1k/1k_ERR1413593_2.fastq.gz"
# aligner_index_dir = "tests/assets/index/bwa"
# chromsizes = "tests/assets/chromsizes/M129.sizes"
# fragment_index = "tests/assets/fragmentIndex/M129_HindIII.bed"


# sample1 = Record(
#     entrypoint=Entrypoint.FASTQ_ALIGN,
#     input_data=Bam(bam=fastq1),#FastqPairedEnd(fastq1=fastq1, fastq2=fastq2),
#     aligner = BWA.MEM2,
#     aligner_index_dir = aligner_index_dir,
#     aligner_index_prefix = "M129",

# )
# print(sample1.model_dump_json())

#     # pairs_label: PairsLabel = lambda data: PairsLabel()
#     # pairs_select: Optional[PairsSelect]
#     # pairs_merge_before_dedup: Optional[PairsMerge]
#     # pairs_dedup: Optional[PairsDedup]
#     # pairs_merge_after_dedup: Optional[PairsMerge]
#     # pairs_hic_bin_coarsen_addnorm: Optional[PairsHiCBinCoarsenAddnorm]
#     # pairs_cool_bin: Optional[PairsCoolBin]
#     # cool_coarsen_addnorm: Optional[CoolCoarsenAddnorm]