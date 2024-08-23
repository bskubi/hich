from collections import Counter
from hich.parse.pairs_file import PairsFile
from numpy import searchsorted
from dataclasses import dataclass, field

@dataclass
class PairsSampler:
    conjuncts: list[str]
    cis_strata: list = None
    event_code: str = ""

    def __post_init__(self):
        cols = ", ".join(self.conjuncts)
        self.event_code = "(" + cols + ")"
        self.event_code = compile(self.event_code, '<string>', 'eval')
        nat_partition = lambda cuts: sorted(list(set(cuts + [float('inf')]))) if cuts else None
        self.cis_strata = nat_partition(self.cis_strata)

    def count(self, pairs_file: PairsFile, lines = float('inf')):
        events = Counter()
        strata = self.cis_strata if self.cis_strata else None

        if strata:
            get_stratum = lambda pair: strata[searchsorted(strata, pair.distance)] if pair.is_cis() else ""
        else:
            get_stratum = lambda pair: None

        for i, pair in enumerate(pairs_file):
            stratum = get_stratum(pair)
            event = eval(self.event_code)
            events[event] += 1
            if i >= lines: break

        return events