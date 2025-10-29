from enum import Enum
from typing import Literal, Optional
from pydantic import BaseModel, Field

from ..config.command import Command

class AnalysisType(str, Enum):
    HICREP = "HICREP",

class BaseAnalysis(BaseModel):
    analysis_type: AnalysisType
    id: str

class HiCRepAnalysis(BaseAnalysis):
    analysis_type: Literal[AnalysisType.HICREP] = AnalysisType.HICREP
    source_ids: list[str]
    resolutions: Optional[list[int]] = None
    h: Optional[list[int]] = None

    command_hich_matrix_hicrep: Command = Field(
            default_factory = lambda data: Command(
            base_command = "hich matrix hicrep",
            opts = {
                "-r": data.get('resolutions', [20_000]),
                "-h": data.get('h', 1),
            },
            flags = [],
            args=["${scc}", "${matrix}"]
        )
    )