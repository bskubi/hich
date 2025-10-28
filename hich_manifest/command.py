from typing import Union, Optional
from typing_extensions import Self
from enum import Enum
from pydantic import BaseModel, Field, FilePath, DirectoryPath, model_validator, field_validator, AliasChoices

class Command(BaseModel):
    base_command: Union[Enum, str]
    opts: dict = Field({})
    flags: list[str] = Field([])
    args: list[str] = Field([])