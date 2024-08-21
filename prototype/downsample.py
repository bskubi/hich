from scipy.optimize import linprog
from hich.parse.pairs_file import PairsFile
from hich.parse.pairs_segment import PairsSegment
import time
import polars as pl
from dataclasses import dataclass, field
from numpy import searchsorted
from collections import defaultdict
from copy import deepcopy
from random import choices
import numpy as np
from abc import ABC, abstractmethod
from types import MethodType 
always = lambda func, vec: all([func(elem) for elem in vec])

@dataclass
class PairsSpace(ABC):
    ignore_code: str = ""

    def __post_init__(self):
        self.ignore_code = compile(self.ignore_code, '<string>', 'eval') if self.ignore_code else None

    def ignore(self, pair): return eval(self.ignore_code)

    @abstractmethod
    def event(self, pair): ...

@dataclass
class TransCisThresholds(PairsSpace):
    cis_thresholds: list = field(default_factory = list)

    def __post_init__(self):
        self.cis_thresholds.sort()
        self.cis_thresholds.append(float('inf'))

    def event(self, pair):
        if pair.is_trans():
            chrom = (pair.chr1, pair.chr2)
            threshold = None
        else:
            chrom = pair.chr1
            threshold = self.cis_threshold(pair.distance)
        return (chrom, threshold)
    
    def cis_threshold(self, distance):
        """Return smallest cis threshold greater than distance"""
        idx = searchsorted(self.cis_thresholds, pair.distance)
        return self.cis_thresholds[idx]
    
    def __str__(self):
        return f"TransCisThresholds.cis_thresholds: {self.cis_thresholds}"

@dataclass
class PairsHistogram:
    space: PairsSpace
    distribution: defaultdict = field(default_factory = lambda: defaultdict(lambda: 0))

    @classmethod
    def center(cls,
               histograms: list['PairsHistogram'],
               center_func: MethodType = np.mean) -> 'PairsHistogram':
        
        support = set()

        for hist in histograms: support.update(hist.support())
        
        PairsHistogram(space = histograms[0].space)
        

    def support(self):
        return [event for event, outcome in self.distribution if outcome > 0]

    def count_pair(self, pair: PairsSegment):
        """Filter pair and increment read bin

        Some pairs will be screened out by the space's ignore_code,
        and the rest will be assigned an event.

        Read bins are defined per chromosome pair (trans) or over a
        uniform set of cis bin thresholds (cis).

        A custom_parser can be a function taking "pair" as an argument and
        returning an event.
        """
        if self.space.ignore(pair): return
        event = self.space.event(pair)
        self.distribution[event] += 1
    
    def events(self): return self.distribution.keys()

    def outcomes(self): return self.distribution.values()

    def total(self): return sum(self.outcomes())

    def set_outcomes(self, new_outcomes):
        """Given list of outcomes, set outcomes in same order as original events"""
        for event, new_outcome in zip(self.events(), new_outcomes):
            self[event] = new_outcome
    
    def is_count(self):
        """Check if count (all outcomes positive integers)"""
        positive_int = lambda outcome: outcome >= 0 and isinstance(outcome, int)
        
        return always(positive_int, self.outcomes())
    
    def is_probdist(self):
        """Check if probability distribution (all outcomes positive, sum to 1)"""
        positive = lambda outcome: outcome >= 0
        outcomes_positive = always(positive, self.outcomes())
        epsilon = .0001

        sum_approx_1 = abs(1.0 - self.total()) < epsilon

        return outcomes_positive and sum_approx_1
               
    def to_probdist(self):
        """Convert to probability distribution (all outcomes positive, sum to 1)"""
        prob = self.copy()
        if not prob.is_probdist():
            for k in prob.distribution:
                prob[k] /= self.total()
        return prob

    def to_count(self, n):
        """Set new read count or fraction of original count, same probabilities"""
        if 0 <= n and n <= 1 and isinstance(n, float):
            new_count = self.total()*n
            rounded = np.round(new_count, 0).as_type(int)
            return self.to_count(rounded)

        count = self.copy()
        outcomes = np.array(list(self.outcomes()))

        # Get ratio of target count to current total count
        downsample = n/count.total()
        
        # Get counts with original probabilities and new total
        fitted = downsample*outcomes

        # Round to integers
        rounded = np.round(fitted, 0).astype(int)

        # Return as new PairsHistogram
        count.set_outcomes(rounded)
        return count

    def downsample_to_probdist(self, prob_dist):
        """Downsample counts to match a probability distribution

        Given a probability distribution, uses linear programming to obtain a
        new number of counts N such that N*prob_dist will minimally downsample
        the original counts histogram, ensuring all counts are positive and
        none are larger than the original count for that event.
        """
        # Ensure the starting counts histogram is in fact positive integers
        assert self.is_count()

        # Cost function (linear programming attempts to minimize cost, so we
        # aim to minimize the negative sum of counts in order to maximize
        # the sum of counts)
        c = [-1]

        # Upper bounds. A_ub is a one-column vector of the probabilities.
        # b_ub is a list of the original counts, serving as the upper bound.
        # Linear programming will ensure that [A_ub]N <= b_ub
        A_ub = np.array([prob_dist]).T
        b_ub = list(self.outcomes())

        # Set lower bound of 0 on N to force positive counts
        bounds = [(0, None)]

        # Find value of N meeting criteria
        result = linprog(c, A_ub=A_ub, b_ub=[b_ub], bounds=bounds)
        assert result.success, "No solution found attempting to match distribution"

        # Get value of N - total number of downsampled events
        downsampled_N = result.x[0]

        # Set downsampled counts and round counts to integers
        count_dist = downsampled_N*np.array(prob_dist)
        rounded_count_dist = np.round(count_dist, 0).astype(int)
        

        # Set outcomes
        pairs_distribution = self.copy()
        pairs_distribution.set_outcomes(rounded_count_dist)

        return pairs_distribution

    def __getitem__(self, key):
        return self.distribution[key]
    
    def __setitem__(self, key, value):
        self.distribution[key] = value
    
    def copy(self): return deepcopy(self)

    def __str__(self):
        s = "\n".join([f"{str(self.space)}",
             f"Distribution: {dict(self.distribution)}"])
        return s



pairs_file = PairsFile("240802_Arima_reseq_dedup.pairs.gz")
cis_thresholds = [10, 20, 50, 100]
space = TransCisThresholds(ignore_code = "not pair.ur()", cis_thresholds = cis_thresholds)

dist = PairsHistogram(space)

for i, pair in enumerate(pairs_file):
    dist.count_pair(pair)
    if i >= 1000:
        break

print(dist)
print(dist.to_probdist())
prob = [.7, .15, .05, .1, 0]
print(dist.downsample_to_probdist(prob))
print(dist.downsample_to_probdist(prob).to_count(50))

"""
0. On input
    Associate samples with downsampling groups
        Default: all in their own target profile group

    Specify a minimum number of reads per contig.
        - Above: follow normal methods
        - Below: either retain all or drop all

1. To form target probdist PairHistograms, convert count PairHistograms to
probabilities and average them, leaving profiles unchanged if they are alone
in the sample group.
    - Can also provide an explicit target for each sample group
    - Specify some sample files and/or contigs to ignore

2. The user can specify:
    2a) (optional) Downsample to target probdist for sample group
    2b) (optional) Downsample to any of:
        - Specific size
        - Specific fraction
        - Minimum sample size within sample group
        - Minimum sample size across sample groups

3. This yields new PairsHistograms for each sample. We then use MultiSampler
to selection sample the pairs files.
"""