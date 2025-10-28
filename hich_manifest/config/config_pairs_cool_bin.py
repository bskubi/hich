from typing import Optional, Annotated
from pydantic import model_validator, field_validator, Field, AfterValidator, FilePath
from .config_process_task_control import ConfigProcessTaskControl
from .command import Command
from .task_control import TaskControl
from .config_contact_matrix import MatrixBinResolutionsValidator, DEFAULT_BIN_RESOLUTIONS


class ConfigPairsCoolBin(ConfigProcessTaskControl):
    genomic_bin_segmentation: Optional[FilePath] = None
    bin_resolutions: Optional[Annotated[list[int], MatrixBinResolutionsValidator]] = DEFAULT_BIN_RESOLUTIONS
    bin_resolution: int = Field(
        default_factory = lambda data: min(data.get('bin_resolutions', [min(DEFAULT_BIN_RESOLUTIONS)]))
    )
    bins: str = Field(
        default_factory = lambda data: (
            data.get('genomic_bin_segmentation').name if data.get('genomic_bin_segmentation') is not None
            else '${chromsizes}:' + str(data.get('bin_resolution'))
        )
    )
    command_cooler_cload_pairs: Command = Field(
        default_factory = lambda data: Command(
            base_command = "cooler cload pairs",
            opts = {
                "-c1": "2",
                "-p1": "3",
                "-c2": "4",
                "-p2": "5",
                "--temp-dir": "."
            },
            args = [data.get('bins'), "${pairs}", "${cool}"]
        )
    )
