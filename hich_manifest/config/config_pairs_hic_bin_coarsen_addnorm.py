from typing import Optional, Annotated
from pydantic import model_validator, Field
from .config_process_task_control import ConfigProcessTaskControl
from .command import Command
from .task_control import TaskControl
from .config_contact_matrix import MatrixBinResolutionsValidator, DEFAULT_BIN_RESOLUTIONS

class ConfigPairsHiCBinCoarsenAddNorm(ConfigProcessTaskControl):
    bin_resolutions: Optional[Annotated[list[int], MatrixBinResolutionsValidator]] = DEFAULT_BIN_RESOLUTIONS
    bin_resolutions_arg: str = Field(default_factory = lambda data: ",".join([str(r) for r in data.get('bin_resolutions', DEFAULT_BIN_RESOLUTIONS)]))

    command_juicer_tools_pre: Command = Field(
        default_factory = lambda data: (
            Command(
                base_command = "juicer_tools pre",
                opts = {
                    "-r": data['bin_resolutions_arg'],
                    "--threads": "${cpus}",
                    "-t": "."
                },
                flags = ["-Xms${xms}g", "-Xmx${xmx}g"],
                args = ["${pairs}", "${hic}", "${chromsizes}"]
            )
        )
    )