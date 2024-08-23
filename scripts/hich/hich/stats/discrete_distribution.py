from collections import Counter
import numpy as np
from numpy.typing import NDArray
from copy import deepcopy
from statistics import mean
from numbers import Number
from functools import reduce
from scipy.optimize import linprog

class DiscreteDistribution(Counter):
    @classmethod
    def mean_mass(cls, distributions: list['DiscreteDistribution']) -> 'DiscreteDistribution':
        combined = reduce(lambda a, b: a+b, distributions)
        return combined/combined.total()

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def to_count(self, n: int) -> 'DiscreteDistribution':
        total = n if isinstance(n, int) else self.total() * n
        new_counts = {event: round(prob*total) for event, prob in self.probabilities().items()}
        return DiscreteDistribution(new_counts)

    def downsample_to_probabilities(self, probdist: 'DiscreteDistribution') -> 'DiscreteDistribution':
        """Given a probability distribution, minimally downsample counts to match it.

        Find the maximum N such that 0 <= N*probdist[event] <= self[event] for all events
        """
        # Cost function (linear programming attempts to minimize cost, so we
        # aim to minimize the negative sum of counts in order to maximize
        # the sum of counts)
        c = [-1]

        # Upper bounds. A_ub is a one-column vector of the probabilities.
        # b_ub is a list of the original counts, serving as the upper bound.
        # Linear programming will ensure that [A_ub]N <= b_ub

        events = self.events()
        A_ub = np.array([probdist.outcomes(events)]).T
        b_ub = list(self.outcomes(events))

        # Set lower bound of 0 on N to force positive counts
        bounds = [(0, None)]

        # Find value of N meeting criteria
        result = linprog(c, A_ub=A_ub, b_ub=[b_ub], bounds=bounds)
        
        assert result.success, "No solution found attempting to match distribution"

        # Get value of N - total number of downsampled events
        N = round(result.x[0])

        return self.to_count(N)

    def probabilities(self) -> NDArray[np.float_]:
        probs = self.copy()
        total = probs.total()
        for event in probs: probs[event] /= total
        return probs

    def __lt__(self, other: 'DiscreteDistribution') -> bool: return self.total() < other.total()

    def __le__(self, other: 'DiscreteDistribution') -> bool: return self.total() <= other.total()

    def __gt__(self,other: 'DiscreteDistribution') -> bool: return self.total() > other.total()

    def __ge__(self, other: 'DiscreteDistribution') -> bool: return self.total() >= other.total()

    def __add__(self, other: 'DiscreteDistribution') -> 'DiscreteDistribution':
        support = set(self.events() + other.events())
        event_outcome = {event: self[event] + other[event] for event in support}
        return DiscreteDistribution(event_outcome)
    
    def __truediv__(self, denom: Number) -> 'DiscreteDistribution':
        result = self.copy()
        for event in self:
            result[event] /= denom
        return result

    def copy(self) -> 'DiscreteDistribution': return deepcopy(self)

    def outcomes(self, events = None) -> list[Number]:
        if not events: return list(self.values())
        else: return [self[event] for event in events]

    def events(self) -> list[object]: return list(self.keys())