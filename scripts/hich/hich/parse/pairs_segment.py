class PairsSegment:
    def __init__(self, **kwargs):
        for k, v in kwargs.items():
            setattr(self, k, v)
    
    def to_string(self):
        return "\t".join(self.__dict__.values())
    
    def __str__(self):
        return self.to_string()