from scipy.optimize import linprog
from hich.parse.pairs_file import PairsFile
from hich.parse.pairs_segment import PairsSegment
from hich.coverage.pairs_space import PairsSpace
from dataclasses import dataclass, field
from numpy import searchsorted
from collections import defaultdict
from copy import deepcopy
from random import choices
import numpy as np
from types import MethodType
from numbers import Number

always = lambda func, vec: all([func(elem) for elem in vec])

@dataclass
class PairsHistogram:
    space: PairsSpace
    distribution: defaultdict = field(default_factory = lambda: defaultdict(lambda: 0))

    @classmethod
    def center(cls,
               histograms: list['PairsHistogram'],
               center_func: MethodType = np.mean,
               ignore_events: list[object] = []) -> 'PairsHistogram':
        """Compute central probability distribution a set of histograms
        
        histograms: a list of probability distributions
        center_func: histograms[event] -> central outcome

        Returns a PairsHistogram with central outcomes according to center_func
        """
        # Skip processing if we have one histogram being set to the mean.
        if center_func in [None, np.mean] and len(histograms) == 1:
            return histograms[0].to_probdist()
        
        support = set()
        for hist in histograms: support.update(hist.support())
        set.difference_update(ignore_events)
        
        central = PairsHistogram(space = histograms[0].space)
        for event in support:
            outcomes = [hist[event] for hist in histograms]
            central[event] = center_func(outcomes)
        return central.to_probdist()

    def support(self) -> list[object]:
        return [event for event, outcome in self.distribution if outcome > 0]

    def count_pair(self, pair: PairsSegment) -> None:
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

    def set_outcomes(self, new_outcomes: list[Number]) -> None:
        """Given list of outcomes, set outcomes in same order as original events"""
        for event, new_outcome in zip(self.events(), new_outcomes):
            self[event] = new_outcome
    
    def is_count(self) -> bool:
        """Check if count (all outcomes positive integers)"""
        positive_int = lambda outcome: outcome >= 0 and isinstance(outcome, int)
        
        return always(positive_int, self.outcomes())
    
    def is_probdist(self) -> bool:
        """Check if probability distribution (all outcomes positive, sum to 1)"""
        positive = lambda outcome: outcome >= 0
        outcomes_positive = always(positive, self.outcomes())
        epsilon = .0001

        sum_approx_1 = abs(1.0 - self.total()) < epsilon

        return outcomes_positive and sum_approx_1
               
    def to_probdist(self) -> 'PairsHistogram':
        """Convert to probability distribution (all outcomes positive, sum to 1)"""
        prob = self.copy()
        if not prob.is_probdist():
            for k in prob.distribution:
                prob[k] /= self.total()
        return prob

    def to_count(self, n: int) -> 'PairsHistogram':
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

    def downsample_to_probdist(self, prob_dist: 'PairsHistogram') -> 'PairsHistogram':
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

    def __getitem__(self, key: object) -> Number:
        return self.distribution[key]
    
    def __setitem__(self, key: object, value: Number) -> None:
        self.distribution[key] = value
    
    def copy(self) -> 'PairsHistogram': return deepcopy(self)

    def __str__(self) -> str:
        s = "\n".join([f"{str(self.space)}",
             f"Distribution: {dict(self.distribution)}"])
        return s