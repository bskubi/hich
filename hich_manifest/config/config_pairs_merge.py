from typing import Union, Optional
from enum import Enum
from pydantic import model_validator, Field

from .config_process import ConfigProcess
from .command import Command


class ConfigPairsMerge(ConfigProcess):
    source_ids: list[str] = Field(min_length=1)
    
    command_pairtools_merge: Command = Command(
        base_command = "pairtools merge",
        opts = {
            "--output": "${pairs_merged}",
            "--nproc-in": "${cpus}",
            "--nproc-out": "${cpus}",
            "--tmpdir": "."
        },
        args = ["${pairs}"]
    )