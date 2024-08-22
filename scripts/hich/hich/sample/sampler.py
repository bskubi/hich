from dataclasses import dataclass, field
from random import random

@dataclass
class Sampler:
    N_records: int
    n_to_sample: int
    t_viewed: int = 0
    m_sampled: int = 0

    def sample(self):
        u = random()
        sample_it = (self.N_records - self.t_viewed) * u < self.n_to_sample - self.m_sampled
        self.t_viewed += 1
        self.m_sampled += sample_it
        return sample_it
    
    def finished(self): return self.m_sampled >= self.n_to_sample
