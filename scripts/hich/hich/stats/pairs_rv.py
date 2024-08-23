from hich.parse.pairs_segment import PairsSegment
from hich.stats.discrete_distribution import DiscreteDistribution
from dataclasses import dataclass
from typing import Callable
from numpy import searchsorted
from polars import DataFrame

def nat_partition(cuts):
    return sorted(list(set(cuts + [float('inf')]))) if cuts else None

@dataclass
class PairsRV:
    conjuncts: list[str] = None
    cis_strata: list = None
    event_code: str = ""

    def __post_init__(self):
        cols = ", ".join(self.conjuncts)
        self.event_code = "(" + cols + ")"
        self.compiled_event_code = compile(self.event_code, '<string>', 'eval')
        
        self.cis_strata = nat_partition(self.cis_strata)

    def get_stratum(self, pair):
        if self.cis_strata:
            return self.cis_strata[searchsorted(self.cis_strata, pair.distance)] if pair.is_cis() else ""
        else:
            return None 

    def event(self, record: PairsSegment):
        strata = self.cis_strata
        stratum = self.get_stratum(record)
        return eval(self.compiled_event_code)
    
    def to_polars(self, distribution: DiscreteDistribution):
        schema = self.conjuncts + ["count"]
        rows = [event + (count,) for event, count in distribution.items()]
        return DataFrame(rows, schema=schema, orient='row')
    
    def from_polars(self, df):
        schema = [col for col in df.columns if col != "count"]
        distribution = DiscreteDistribution()
        for row in df.iter_rows():
            event = row[:-1]
            count = row[-1]
            distribution[event] = count
        return distribution

    def __getstate__(self):
        # Return a dictionary excluding the compiled code
        state = self.__dict__.copy()
        del state['compiled_event_code']
        return state

    def __setstate__(self, state):
        # Restore the object's state and recompile the code
        self.__dict__.update(state)
        self.compiled_event_code = compile(self.event_code, '<string>', 'eval')

