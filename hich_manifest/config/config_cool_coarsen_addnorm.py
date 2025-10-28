from typing import Optional, Annotated
from pydantic import model_validator, Field
from .config_process_task_control import ConfigProcessTaskControl
from .command import Command
from .config_contact_matrix import MatrixBinResolutionsValidator, DEFAULT_BIN_RESOLUTIONS

class ConfigCoolCoarsenAddNorm(ConfigProcessTaskControl):
    bin_resolutions: Optional[Annotated[list[int], MatrixBinResolutionsValidator]] = DEFAULT_BIN_RESOLUTIONS
    bin_resolutions_arg: str = Field(default_factory = lambda data: ",".join([str(r) for r in data.get('bin_resolutions', DEFAULT_BIN_RESOLUTIONS)]))

    command_cooler_zoomify: Command = Field(
        default_factory = lambda data: (
            Command(
                base_command = "cooler zoomify",
                opts = {
                    "--nproc": "${cpus}",
                    "--balance-args": "--nproc ${cpus} --max-iters 2000",
                    "--resolutions": data['bin_resolutions_arg'],
                    "--out": "${mcool}"
                },
                flags = ["--balance"],
                args = ["${mcool}"]
            )
        )
    )
