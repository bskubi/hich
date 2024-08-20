@dataclass
class StratumSampler:
    
    N_records: int
    n_to_sample: int
    t_viewed: int = 0
    m_sampled: int = 0
    U_random: list[float] = field(default_factory=list)
    samples: list[bool] = field(default_factory=list)

    def batch(self):
        """Add a number of sample values equal to the length of self.U_random"""
        for u in self.U_random:
            if (self.N_records - self.t_viewed) * u < self.n_to_sample - self.m_sampled:
                self.m_sampled += 1
                self.samples.append(True)
            else:
                self.samples.append(False)
            self.t_viewed += 1
    
    def sample(self, n):
        """Pop and return up to the first n items from self.samples"""
        popped = self.samples[:n]
        self.samples = self.samples[n:]
        return popped