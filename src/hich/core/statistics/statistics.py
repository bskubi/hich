import polars as pl
from abc import ABC, abstractmethod
from dataclasses import dataclass

@dataclass
class Statistic:
    df: pl.DataFrame

class StatisticType:
    @property
    @abstractmethod
    def accumulator(self) -> type["StatisticAccumulator"]: ...

class StatisticAccumulator(ABC):
    @abstractmethod
    def update(self, obs_df: pl.DataFrame): ...
    
    @property
    @abstractmethod
    def statistic(self) -> Statistic: ...