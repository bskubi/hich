from dataclasses import dataclass, field
from hich.parse.pairs_segment import PairsSegment
from abc import ABC, abstractmethod

@dataclass
class PairsSpace(ABC):
    ignore_code: str = ""

    def __post_init__(self):
        self.ignore_code = compile(self.ignore_code, '<string>', 'eval') if self.ignore_code else None

    def ignore(self, pair: PairsSegment) -> bool: return eval(self.ignore_code)

    @abstractmethod
    def event(self, pair: PairsSegment) -> object: ...
