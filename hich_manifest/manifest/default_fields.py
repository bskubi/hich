from typing import Literal, Annotated, Union, Optional
from pydantic import BaseModel, Field, FilePath, DirectoryPath

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

DEFAULT_BIN_RESOLUTIONS_FIELD = Field(DEFAULT_BIN_RESOLUTIONS, description="Contact matrix binning resolutions (bp)")

DEFAULT_CONFIG_BIN_PARSE_PAIRS_FIELD = Field(
    default_factory = lambda data: ConfigBamParsePairs(chromsizes=data.get('chromsizes'))
)

DEFAULT_CONFIG_PAIRS_LABEL_FIELD = Field(
    default_factory = lambda data: (
        ConfigPairsLabel(task=TaskControl.RUN if data.get('fragment_index') else TaskControl.SKIP)
    )
)

DEFAULT_CONFIG_PAIRS_SELECT_FIELD = Field(
    default_factory = lambda data: (
        ConfigPairsSelect( 
            pairs_filters = PairsFilters(
                not_same_frag=data.get("fragment_index") is not None
            )
        )
    )
)

DEFAULT_CONFIG_PAIRS_MERGE_FIELD = Field(
    default_factory = lambda data: (
        ConfigPairsMerge()
    )
)

DEFAULT_CONFIG_PAIRS_DEDUP_FIELD = ConfigPairsDedup()

DEFAULT_CONFIG_PAIRS_HIC_BIN_COARSEN_ADDNORM_FIELD = Field(
    default_factory = lambda data: ConfigPairsHiCBinCoarsenAddNorm(bin_resolutions = data.get('bin_resolutions'))
)

DEFAULT_CONFIG_PAIRS_COOL_BIN_FIELD = Field(
    default_factory = lambda data: ConfigPairsCoolBin(bin_resolutions = data.get('bin_resolutions'))
)

DEFAULT_CONFIG_COOL_COARSEN_ADDNORM_FIELD = Field(
    default_factory = lambda data: ConfigCoolCoarsenAddNorm(bin_resolutions = data.get('bin_resolutions'))
)