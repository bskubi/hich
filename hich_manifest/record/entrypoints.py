from typing import Union, Optional, Literal
from typing_extensions import Self
from enum import Enum
from pydantic import BaseModel, Field, FilePath, DirectoryPath, model_validator, field_validator, AliasChoices

class Entrypoint(str, Enum):
    FASTQ_ALIGN = "FASTQ_ALIGN"
    BAM_PARSE_PAIRS = "BAM_PARSE_PAIRS"
    PAIRS_LABEL = "PAIRS_LABEL"
    PAIRS_SELECT = "PAIRS_SELECT"
    PAIRS_MERGE_BEFORE_DEDUP = "PAIRS_MERGE_BEFORE_DEDUP"
    PAIRS_DEDUP = "PAIRS_DEDUP"
    PAIRS_MERGE_AFTER_DEDUP = "PAIRS_MERGE_AFTER_DEDUP"
    PAIRS_BIN_COARSEN_ADDNORM = "PAIRS_BIN_COARSEN_ADDNORM"
    COOL_COARSEN_ADDNORM = "COOL_COARSEN_ADDNORM"


