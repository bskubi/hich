class MultiSampler:
    def __init__(self, df, id_cols, N_records_col, n_to_sample_col, seed = 42):
        self.setup(df, id_cols, N_records_col, n_to_sample_col, seed)

    def setup(self, df, id_cols, N_records_col, n_to_sample_col, seed = 42):
        sel = collapse([id_cols, N_records_col, n_to_sample_col])
        df = df.select(sel)

        N_col = len(id_cols)
        n_col = N_col + 1

        self.id_cols = id_cols
        self.rng = default_rng(seed=seed)
        self.samplers = {}
        for row in df.iter_rows():
            sampler_id = row[:len(id_cols)]
            N = row[N_col]
            n = row[n_col]
            self.samplers[sampler_id] = StratumSampler(N, n)
    
    def sample(self, df):
        groups = df.with_row_index() \
                   .group_by(self.id_cols)
        sampled = []
        for name, data in groups:
            n = len(data)
            sample = self.batch(name, n)
            sampled.append(data.filter(sample))
        
        index_col = df.columns[0]
        sampled = pl.concat(sampled).sort(by=index_col)
        return sampled
    
    def batch(self, sampler_id, n):
        self.samplers[sampler_id].U_random += self.rng.random(n).tolist()
        self.samplers[sampler_id].batch()
        return self.samplers[sampler_id].sample(n)