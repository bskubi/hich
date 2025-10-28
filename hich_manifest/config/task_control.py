from enum import Enum

class TaskControl(str, Enum):
    RUN = "RUN"
    SKIP = "SKIP"
