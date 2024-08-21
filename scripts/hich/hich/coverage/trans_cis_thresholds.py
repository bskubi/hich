from dataclasses import dataclass, field
from hich.coverage.pairs_space import PairsSpace
from numpy import searchsorted
from hich.parse.pairs_segment import PairsSegment
from numbers import Number

@dataclass
class TransCisThresholds(PairsSpace):
    cis_thresholds: list = field(default_factory = list)

    def __post_init__(self):
        self.cis_thresholds.sort()
        self.cis_thresholds.append(float('inf'))

    def event(self, pair: PairsSegment) -> object:
        if pair.is_trans():
            chrom = (pair.chr1, pair.chr2)
            threshold = None
        else:
            chrom = pair.chr1
            threshold = self.cis_threshold(pair.distance)
        return (chrom, threshold)
    
    def cis_threshold(self, distance: Number) -> Number:
        """Return smallest cis threshold greater than distance"""
        idx = searchsorted(self.cis_thresholds, distance)
        return self.cis_thresholds[idx]
    
    def __str__(self):
        return f"TransCisThresholds.cis_thresholds: {self.cis_thresholds}"