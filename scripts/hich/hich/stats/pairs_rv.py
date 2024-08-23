from hich.parse.pairs_segment import PairsSegment

class PairsRV:
    conjuncts: list[str]
    cis_strata: list = None
    event_code: str = ""

    def __post_init__(self):
        cols = ", ".join(self.conjuncts)
        self.event_code = "(" + cols + ")"
        self.event_code = compile(self.event_code, '<string>', 'eval')
        nat_partition = lambda cuts: sorted(list(set(cuts + [float('inf')]))) if cuts else None
        self.cis_strata = nat_partition(self.cis_strata)
        
        if self.cis_strata:
            self.get_stratum = lambda pair: strata[searchsorted(strata, pair.distance)] if pair.is_cis() else ""
        else:
            self.get_stratum = lambda pair: None

    def event(self, record: PairsSegment):
        strata = self.cis_strata
        stratum = get_stratum(record)
        return eval(self.event_code)