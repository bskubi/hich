from typing import Optional
from pydantic import model_validator
from .config_process_task_control import ConfigProcessTaskControl
from .command import Command
from .task_control import TaskControl


class ConfigPairsDedup(ConfigProcessTaskControl):
    command_pairtools_dedup: Command = Command(
            base_command = "pairtools dedup",
            opts = {
                "--output": "${pairs_deduped}",
                "--nproc-in": "${cpus}",
                "--nproc-out": "${cpus}"
            },
            args = ["${pairs}"]
        )