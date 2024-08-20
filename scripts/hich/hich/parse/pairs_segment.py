class PairsSegment:
    reserved = {"readID": str, "chr1": str, "pos1": int, "chr2": str, "pos2": int, "strand1": str, "strand2": str}
    required = {"chr1": str, "pos1": int, "chr2": str, "pos2": int}
    alt = {"chrom1":"chr1", "chrom2":"chr2"}
    extra = {"distance": lambda pair: setattr(pair, "distance", abs(pair.pos1 - pair.pos2) if pair.is_cis() else None)}

    def __init__(self, **kwargs):
        self.__dict__ = {k: None for k in PairsSegment.reserved}
        for k, v in kwargs.items():
            if k in PairsSegment.reserved:
                cast = PairsSegment.reserved[k]
                v = cast(v)
            setattr(self, k,v)
        for alt, main in PairsSegment.alt.items():
            print(main not in self.__dict__, alt in self.__dict)

            if alt not in self.__dict__ and main in self.__dict__:
                setattr(self, alt, self.__dict__[main])
            elif main not in self.__dict__ and alt in self.__dict__:
                setattr(self, main, self.__dict__[alt])
        for func in PairsSegment.extra.values():
            func(self)

    def to_dict(self, columns = None):
        if columns:
            return {c:self.__dict__[c] for c in columns if c in self.__dict__}
        else:
            return {k:v for k, v in self.__dict__.items() if k not in PairsSegment.alt}

    def to_string(self, columns = None):
        return "\t".join(str(v) for v in self.to_dict(columns).values())

    def is_cis(self): return self.chr1 == self.chr2
    
    def is_trans(self): return self.chr1 != self.chr2

    def intrachrom(self): return self.is_cis()

    def interchrom(self): return self.is_trans()

    def ur(self): return self.pair_type in ["UR", "RU", "UR"]