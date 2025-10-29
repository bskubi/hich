from typing import Literal, Annotated, Union, Optional
from pydantic import BaseModel, Field, FilePath, DirectoryPath, RootModel

from .entrypoints import Entrypoint
from ..config.task_control import TaskControl
from ..config.config_fastq_align import ConfigFastqAlign
from ..config.config_bam_parse_pairs import ConfigBamParsePairs
from ..config.config_pairs_label import ConfigPairsLabel
from ..config.config_pairs_select import ConfigPairsSelect, PairsFilters
from ..config.config_pairs_merge import ConfigPairsMerge
from ..config.config_pairs_dedup import ConfigPairsDedup
from ..config.config_pairs_hic_bin_coarsen_addnorm import ConfigPairsHiCBinCoarsenAddNorm
from ..config.config_pairs_cool_bin import ConfigPairsCoolBin
from ..config.config_cool_coarsen_addnorm import ConfigCoolCoarsenAddNorm
from ..config.config_contact_matrix import MatrixBinResolutionsValidator, DEFAULT_BIN_RESOLUTIONS

from .default_fields import *

class BaseRecord(BaseModel):
    id: Optional[str] = None
    entrypoint: Entrypoint

class FastqAlignRecord(BaseRecord):
    entrypoint: Literal[Entrypoint.FASTQ_ALIGN]
    fastq1: FilePath
    fastq2: Optional[FilePath] = None
    aligner_index_dir: DirectoryPath
    chromsizes: FilePath
    fragment_index: Optional[FilePath] = None
    chrom_subset: Optional[FilePath] = None
    bin_resolutions: Optional[Annotated[list[int], MatrixBinResolutionsValidator]] = DEFAULT_BIN_RESOLUTIONS_FIELD

    config_fastq_align: ConfigFastqAlign
    config_bam_parse_pairs: ConfigBamParsePairs = DEFAULT_CONFIG_BIN_PARSE_PAIRS_FIELD
    config_pairs_label: ConfigPairsLabel = DEFAULT_CONFIG_PAIRS_LABEL_FIELD
    config_pairs_select: ConfigPairsSelect = DEFAULT_CONFIG_PAIRS_SELECT_FIELD
    config_pairs_dedup: ConfigPairsDedup = DEFAULT_CONFIG_PAIRS_DEDUP_FIELD
    config_pairs_hic_bin_coarsen_addnorm: ConfigPairsHiCBinCoarsenAddNorm = DEFAULT_CONFIG_PAIRS_HIC_BIN_COARSEN_ADDNORM_FIELD
    config_pairs_cool_bin: ConfigPairsCoolBin = DEFAULT_CONFIG_PAIRS_COOL_BIN_FIELD
    config_cool_coarsen_addnorm: ConfigCoolCoarsenAddNorm = DEFAULT_CONFIG_COOL_COARSEN_ADDNORM_FIELD


class BamParsePairsRecord(BaseRecord):
    entrypoint: Literal[Entrypoint.BAM_PARSE_PAIRS]
    bam: FilePath
    chromsizes: FilePath
    fragment_index: Optional[FilePath] = None
    chrom_subset: Optional[FilePath] = None
    bin_resolutions: Optional[Annotated[list[int], MatrixBinResolutionsValidator]] = DEFAULT_BIN_RESOLUTIONS_FIELD

    config_bam_parse_pairs: ConfigBamParsePairs = DEFAULT_CONFIG_BIN_PARSE_PAIRS_FIELD
    config_pairs_label: ConfigPairsLabel = DEFAULT_CONFIG_PAIRS_LABEL_FIELD
    config_pairs_select: ConfigPairsSelect = DEFAULT_CONFIG_PAIRS_SELECT_FIELD
    config_pairs_dedup: ConfigPairsDedup = DEFAULT_CONFIG_PAIRS_DEDUP_FIELD
    config_pairs_hic_bin_coarsen_addnorm: ConfigPairsHiCBinCoarsenAddNorm = DEFAULT_CONFIG_PAIRS_HIC_BIN_COARSEN_ADDNORM_FIELD
    config_pairs_cool_bin: ConfigPairsCoolBin = DEFAULT_CONFIG_PAIRS_COOL_BIN_FIELD
    config_cool_coarsen_addnorm: ConfigCoolCoarsenAddNorm = DEFAULT_CONFIG_COOL_COARSEN_ADDNORM_FIELD

class PairsLabelRecord(BaseRecord):
    entrypoint: Literal[Entrypoint.PAIRS_LABEL]
    pairs: FilePath

    chromsizes: FilePath
    fragment_index: Optional[FilePath] = None
    chrom_subset: Optional[FilePath] = None
    bin_resolutions: Optional[Annotated[list[int], MatrixBinResolutionsValidator]] = DEFAULT_BIN_RESOLUTIONS_FIELD

    config_pairs_label: ConfigPairsLabel = DEFAULT_CONFIG_PAIRS_LABEL_FIELD
    config_pairs_select: ConfigPairsSelect = DEFAULT_CONFIG_PAIRS_SELECT_FIELD
    config_pairs_dedup: ConfigPairsDedup = DEFAULT_CONFIG_PAIRS_DEDUP_FIELD
    config_pairs_hic_bin_coarsen_addnorm: ConfigPairsHiCBinCoarsenAddNorm = DEFAULT_CONFIG_PAIRS_HIC_BIN_COARSEN_ADDNORM_FIELD
    config_pairs_cool_bin: ConfigPairsCoolBin = DEFAULT_CONFIG_PAIRS_COOL_BIN_FIELD
    config_cool_coarsen_addnorm: ConfigCoolCoarsenAddNorm = DEFAULT_CONFIG_COOL_COARSEN_ADDNORM_FIELD

class PairsSelectRecord(BaseRecord):
    entrypoint: Literal[Entrypoint.PAIRS_SELECT]
    pairs: FilePath

    chromsizes: FilePath
    chrom_subset: Optional[FilePath] = None
    bin_resolutions: Optional[Annotated[list[int], MatrixBinResolutionsValidator]] = DEFAULT_BIN_RESOLUTIONS_FIELD

    config_pairs_select: ConfigPairsSelect = DEFAULT_CONFIG_PAIRS_SELECT_FIELD
    config_pairs_dedup: ConfigPairsDedup = DEFAULT_CONFIG_PAIRS_DEDUP_FIELD
    config_pairs_hic_bin_coarsen_addnorm: ConfigPairsHiCBinCoarsenAddNorm = DEFAULT_CONFIG_PAIRS_HIC_BIN_COARSEN_ADDNORM_FIELD
    config_pairs_cool_bin: ConfigPairsCoolBin = DEFAULT_CONFIG_PAIRS_COOL_BIN_FIELD
    config_cool_coarsen_addnorm: ConfigCoolCoarsenAddNorm = DEFAULT_CONFIG_COOL_COARSEN_ADDNORM_FIELD

class PairsMergeBeforeDedupSourceRecord(BaseRecord):
    entrypoint: Literal[Entrypoint.PAIRS_MERGE_BEFORE_DEDUP_SOURCE]
    pairs: FilePath

    chromsizes: FilePath
    bin_resolutions: Optional[Annotated[list[int], MatrixBinResolutionsValidator]] = DEFAULT_BIN_RESOLUTIONS_FIELD

    config_pairs_dedup: ConfigPairsDedup = DEFAULT_CONFIG_PAIRS_DEDUP_FIELD
    config_pairs_hic_bin_coarsen_addnorm: ConfigPairsHiCBinCoarsenAddNorm = DEFAULT_CONFIG_PAIRS_HIC_BIN_COARSEN_ADDNORM_FIELD
    config_pairs_cool_bin: ConfigPairsCoolBin = DEFAULT_CONFIG_PAIRS_COOL_BIN_FIELD
    config_cool_coarsen_addnorm: ConfigCoolCoarsenAddNorm = DEFAULT_CONFIG_COOL_COARSEN_ADDNORM_FIELD

class PairsMergeBeforeDedupTargetRecord(BaseRecord):
    entrypoint: Literal[Entrypoint.PAIRS_MERGE_BEFORE_DEDUP_TARGET]

    chromsizes: FilePath
    bin_resolutions: Optional[Annotated[list[int], MatrixBinResolutionsValidator]] = DEFAULT_BIN_RESOLUTIONS_FIELD

    config_pairs_merge: ConfigPairsMerge = DEFAULT_CONFIG_PAIRS_MERGE_FIELD
    config_pairs_dedup: ConfigPairsDedup = DEFAULT_CONFIG_PAIRS_DEDUP_FIELD
    config_pairs_hic_bin_coarsen_addnorm: ConfigPairsHiCBinCoarsenAddNorm = DEFAULT_CONFIG_PAIRS_HIC_BIN_COARSEN_ADDNORM_FIELD
    config_pairs_cool_bin: ConfigPairsCoolBin = DEFAULT_CONFIG_PAIRS_COOL_BIN_FIELD
    config_cool_coarsen_addnorm: ConfigCoolCoarsenAddNorm = DEFAULT_CONFIG_COOL_COARSEN_ADDNORM_FIELD

class PairsDedupRecord(BaseRecord):
    entrypoint: Literal[Entrypoint.PAIRS_DEDUP]
    pairs: FilePath

    chromsizes: FilePath
    bin_resolutions: Optional[Annotated[list[int], MatrixBinResolutionsValidator]] = DEFAULT_BIN_RESOLUTIONS_FIELD

    config_pairs_dedup: ConfigPairsDedup = DEFAULT_CONFIG_PAIRS_DEDUP_FIELD
    config_pairs_hic_bin_coarsen_addnorm: ConfigPairsHiCBinCoarsenAddNorm = DEFAULT_CONFIG_PAIRS_HIC_BIN_COARSEN_ADDNORM_FIELD
    config_pairs_cool_bin: ConfigPairsCoolBin = DEFAULT_CONFIG_PAIRS_COOL_BIN_FIELD
    config_cool_coarsen_addnorm: ConfigCoolCoarsenAddNorm = DEFAULT_CONFIG_COOL_COARSEN_ADDNORM_FIELD

class PairsMergeAfterDedupSourceRecord(BaseRecord):
    entrypoint: Literal[Entrypoint.PAIRS_MERGE_AFTER_DEDUP_SOURCE]
    pairs: FilePath

    chromsizes: FilePath
    bin_resolutions: Optional[Annotated[list[int], MatrixBinResolutionsValidator]] = DEFAULT_BIN_RESOLUTIONS_FIELD

    config_pairs_hic_bin_coarsen_addnorm: ConfigPairsHiCBinCoarsenAddNorm = DEFAULT_CONFIG_PAIRS_HIC_BIN_COARSEN_ADDNORM_FIELD
    config_pairs_cool_bin: ConfigPairsCoolBin = DEFAULT_CONFIG_PAIRS_COOL_BIN_FIELD
    config_cool_coarsen_addnorm: ConfigCoolCoarsenAddNorm = DEFAULT_CONFIG_COOL_COARSEN_ADDNORM_FIELD

class PairsMergeAfterDedupTargetRecord(BaseRecord):
    entrypoint: Literal[Entrypoint.PAIRS_MERGE_AFTER_DEDUP_TARGET]

    chromsizes: FilePath
    bin_resolutions: Optional[Annotated[list[int], MatrixBinResolutionsValidator]] = DEFAULT_BIN_RESOLUTIONS_FIELD

    config_pairs_merge: ConfigPairsMerge = DEFAULT_CONFIG_PAIRS_MERGE_FIELD
    config_pairs_hic_bin_coarsen_addnorm: ConfigPairsHiCBinCoarsenAddNorm = DEFAULT_CONFIG_PAIRS_HIC_BIN_COARSEN_ADDNORM_FIELD
    config_pairs_cool_bin: ConfigPairsCoolBin = DEFAULT_CONFIG_PAIRS_COOL_BIN_FIELD
    config_cool_coarsen_addnorm: ConfigCoolCoarsenAddNorm = DEFAULT_CONFIG_COOL_COARSEN_ADDNORM_FIELD

class PairsBinCoarsenAddNormRecord(BaseRecord):
    entrypoint: Literal[Entrypoint.PAIRS_BIN_COARSEN_ADDNORM]
    pairs: FilePath

    chromsizes: FilePath
    bin_resolutions: Optional[Annotated[list[int], MatrixBinResolutionsValidator]] = DEFAULT_BIN_RESOLUTIONS_FIELD

    config_pairs_hic_bin_coarsen_addnorm: ConfigPairsHiCBinCoarsenAddNorm = DEFAULT_CONFIG_PAIRS_HIC_BIN_COARSEN_ADDNORM_FIELD
    config_pairs_cool_bin: ConfigPairsCoolBin = DEFAULT_CONFIG_PAIRS_COOL_BIN_FIELD
    config_cool_coarsen_addnorm: ConfigCoolCoarsenAddNorm = DEFAULT_CONFIG_COOL_COARSEN_ADDNORM_FIELD

class CoolCoarsenAddNormRecord(BaseRecord):
    entrypoint: Literal[Entrypoint.COOL_COARSEN_ADDNORM]
    cool: FilePath

    chromsizes: FilePath
    bin_resolutions: Optional[Annotated[list[int], MatrixBinResolutionsValidator]] = DEFAULT_BIN_RESOLUTIONS_FIELD

    config_pairs_hic_bin_coarsen_addnorm: ConfigPairsHiCBinCoarsenAddNorm = DEFAULT_CONFIG_PAIRS_HIC_BIN_COARSEN_ADDNORM_FIELD
    config_pairs_cool_bin: ConfigPairsCoolBin = DEFAULT_CONFIG_PAIRS_COOL_BIN_FIELD
    config_cool_coarsen_addnorm: ConfigCoolCoarsenAddNorm = DEFAULT_CONFIG_COOL_COARSEN_ADDNORM_FIELD

RecordUnion = Annotated[
    Union[
        FastqAlignRecord,
        BamParsePairsRecord,
        PairsLabelRecord,
        PairsSelectRecord,
        PairsMergeBeforeDedupSourceRecord,
        PairsMergeBeforeDedupTargetRecord,
        PairsDedupRecord,
        PairsMergeAfterDedupSourceRecord,
        PairsMergeAfterDedupTargetRecord,
        PairsBinCoarsenAddNormRecord,
        CoolCoarsenAddNormRecord
    ],
    Field(discriminator="entrypoint")
]

class Record(RootModel):
    root: RecordUnion
    def __init__(self, **data):
        super().__init__(root=data)
    
    def __getattr__(self, item):
        return getattr(self.root, item)