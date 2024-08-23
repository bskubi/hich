from hich.stats.discrete_distribution import DiscreteDistribution
from dataclasses import dataclass, field
from random import random

@dataclass
class SelectionSampler:
    total: DiscreteDistribution
    target: DiscreteDistribution
    viewed: DiscreteDistribution
    kept: DiscreteDistribution

    def count(self, event):
        self.total[event] += 1

    def sample(self, event):
        u_random = random()
        unseen = (self.total[event] - self.viewed[event])
        undersampled = self.target[event] - self.kept[event]
        keep = unseen * u_random < undersampled
        self.kept[event] += keep
        self.viewed[event] += 1
        return keep

    
