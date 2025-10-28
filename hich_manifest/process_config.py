from pydantic import BaseModel

class ProcessConfig(BaseModel):
    skip: bool = False