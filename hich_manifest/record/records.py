from typing import Literal, Annotated, Union, Optional
from pydantic import BaseModel, Field, FilePath, DirectoryPath

from .entrypoints import Entrypoint
from ..config.task_control import TaskControl
from ..config.config_fastq_align import ConfigFastqAlign
from ..config.config_bam_parse_pairs import ConfigBamParsePairs
from ..config.config_pairs_label import ConfigPairsLabel
from ..config.config_pairs_select import ConfigPairsSelect, PairsFilters
from ..config.config_pairs_dedup import ConfigPairsDedup
from ..config.config_pairs_hic_bin_coarsen_addnorm import ConfigPairsHiCBinCoarsenAddNorm
from ..config.config_pairs_cool_bin import ConfigPairsCoolBin
from ..config.config_cool_coarsen_addnorm import ConfigCoolCoarsenAddNorm
from ..config.config_contact_matrix import MatrixBinResolutionsValidator, DEFAULT_BIN_RESOLUTIONS

class BaseRecord(BaseModel):
    id: str
    entrypoint: Entrypoint

class FastqEntryRecord(BaseRecord):
    bin_resolutions: Optional[Annotated[list[int], MatrixBinResolutionsValidator]] = Field(DEFAULT_BIN_RESOLUTIONS, description="Contact matrix binning resolutions (bp)")
    entrypoint: Literal[Entrypoint.FASTQ_ALIGN]
    fastq1: FilePath
    fastq2: Optional[FilePath] = None
    aligner_index_dir: DirectoryPath
    chromsizes: FilePath
    fragment_index: Optional[FilePath] = None
    chrom_subset: Optional[FilePath] = None

    config_fastq_align: ConfigFastqAlign
    config_bam_parse_pairs: ConfigBamParsePairs = Field(
        default_factory = lambda data: ConfigBamParsePairs(chromsizes=data.get('chromsizes'))
    )
    config_pairs_label: ConfigPairsLabel = Field(
        default_factory = lambda data: (
            ConfigPairsLabel(task=TaskControl.RUN if data.get('fragment_index') else TaskControl.SKIP)
        )
    )
    config_pairs_select: ConfigPairsSelect = Field(
        default_factory = lambda data: (
            ConfigPairsSelect( 
                pairs_filters = PairsFilters(
                    not_same_frag=data.get("fragment_index") is not None
                )
            )
        )
    )
    config_pairs_dedup: ConfigPairsDedup = ConfigPairsDedup()
    config_pairs_hic_bin_coarsen_addnorm: ConfigPairsHiCBinCoarsenAddNorm = Field(
        default_factory = lambda data: ConfigPairsHiCBinCoarsenAddNorm(bin_resolutions = data.get('bin_resolutions'))
    )
    config_pairs_cool_bin: ConfigPairsCoolBin = Field(
        default_factory = lambda data: ConfigPairsCoolBin(bin_resolutions = data.get('bin_resolutions'))
    )
    config_cool_coarsen_addnorm: ConfigCoolCoarsenAddNorm = Field(
        default_factory = lambda data: ConfigCoolCoarsenAddNorm(bin_resolutions = data.get('bin_resolutions'))
    )


class BamParsePairsRecord(BaseRecord):
    entrypoint: Literal[Entrypoint.BAM_PARSE_PAIRS]
    bam: FilePath

class PairsLabel(BaseRecord):
    entrypoint: Literal[Entrypoint.PAIRS_LABEL]
    pairs: FilePath

class PairsSelect(BaseRecord):
    entrypoint: Literal[Entrypoint.PAIRS_SELECT]
    pairs: FilePath

class PairsMergeBeforeDedupSource(BaseRecord):
    entrypoint: Literal[Entrypoint.PAIRS_MERGE_BEFORE_DEDUP]
    pairs: FilePath

class PairsMergeBeforeDedupTarget(BaseRecord):
    entrypoint: Literal[Entrypoint.PAIRS_MERGE_BEFORE_DEDUP]
    source_ids: list[str]

class PairsDedup(BaseRecord):
    entrypoint: Literal[Entrypoint.PAIRS_DEDUP]
    pairs: FilePath

class PairsMergeAfterDedupSource(BaseRecord):
    entrypoint: Literal[Entrypoint.PAIRS_MERGE_AFTER_DEDUP]
    pairs: FilePath

class PairsMergeAfterDedupTarget(BaseRecord):
    entrypoint: Literal[Entrypoint.PAIRS_MERGE_AFTER_DEDUP]
    source_ids: list[str]

class PairsBinCoarsenAddNorm(BaseRecord):
    entrypoint: Literal[Entrypoint.PAIRS_BIN_COARSEN_ADDNORM]
    pairs: FilePath

class CoolCoarsenAddNorm(BaseRecord):
    entrypoint: Literal[Entrypoint.COOL_COARSEN_ADDNORM]
    cool: FilePath

Record = Annotated[
    Union[FastqEntryRecord],
    Field(discriminator="entrypoint")
]