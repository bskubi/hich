from enum import Enum
from typing import Optional
from pydantic import model_validator
from .config_process_task_control import ConfigProcessTaskControl
from .command import Command
from .task_control import TaskControl

class Label(str, Enum):
    RESTRICTION_FRAGMENT = "RESTRICTION_FRAGMENT"


class ConfigPairsLabel(ConfigProcessTaskControl):
    command_hich_pairs_map_ends: Optional[Command] = None

    @model_validator(mode='before')
    @classmethod
    def set_defaults(cls, data: dict) -> dict:
        if not data.get('command_hich_pairs_map_ends') and not data.get('task') == TaskControl.SKIP:
            data['command_hich_pairs_map_ends'] = Command(
                base_command="hich pairs map-ends",
                opts={
                    "--idx1": "rfrag1",
                    "--start1": "rfrag_start1",
                    "--end1": "rfrag_end1",
                    "--idx2": "rfrag2",
                    "--start2": "rfrag_start2",
                    "--end2": "rfrag_end2"
                },
                args=["${fragment_index}", "${pairs}", "${pairs_labeled}"]
            )
        return data