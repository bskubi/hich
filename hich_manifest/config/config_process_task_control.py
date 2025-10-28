from .config_process import ConfigProcess
from .command import Command
from .task_control import TaskControl

from .config_process import ConfigProcess

class ConfigProcessTaskControl(ConfigProcess):
    task: TaskControl = TaskControl.RUN