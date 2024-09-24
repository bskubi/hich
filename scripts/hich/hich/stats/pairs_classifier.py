from dataclasses import dataclass
from hich.parse.pairs_segment import PairsSegment
from hich.stats.discrete_distribution import DiscreteDistribution
from polars import DataFrame
from typing import Callable
import time
import bisect

def nat_partition(cuts):
    return sorted(list(set(cuts + [float('inf')]))) if cuts else None

@dataclass
class PairsClassifier:
    """Classifies PairsSegment objects into outcomes.

    The PairsClassifier defines a way of mapping PairsSegment objects into
    tuples, representing outcomes. An immersive example is:

    PairsClassifier(['record.chr1', 'record.chr2', 'stratum'], [1000, 10000])

    This PairsClassifier will categorize PairsSegment objects ('record') according
    to chr1, chr2, and strata (a binning of the distance between pos1 and pos2).
    If the distance is larger than the highest strata, 'inf' will be returned.

    This works by generating Python code expressing a tuple that combines
    the conjuncts 'record.chr1', 'record.chr2', and 'stratum'. This is compiled
    and used to obtain the relevant data as a tuple. The resulting tuple is
    the return value of the 'classify' function.
    """
    conjuncts: list[str] = None
    cis_strata: list = None
    classification_code: str = ""

    def __post_init__(self):
        """Creates and compiles Python classification code that converts conjuncts into a tuple"""

        # Create a comma-separated list of the conjuncts and add parenthesis
        prefix = lambda conjunct: not (conjunct.strip().startswith("record.") or conjunct in ["strata", "stratum"])
        ensure_format = lambda conjunct: f"record.{conjunct.strip()}" if prefix(conjunct) else conjunct
        formatted_conjuncts = [ensure_format(conjunct) for conjunct in self.conjuncts]
        tuple_center = ", ".join(formatted_conjuncts) if formatted_conjuncts else None
        self.classification_code = "(" + tuple_center + ",)" if tuple_center else ""

        # Compile to Python bytecode
        self.compiled_classification_code = compile(self.classification_code, '<string>', 'eval') if self.classification_code else None

        # If strata is given, add 'inf' to it to clearly define the largest stratum
        self.cis_strata = nat_partition(self.cis_strata)

    def get_stratum(self, pair: PairsSegment):
        """Categorize PairsSegment into stratum
        
            If strata are defined, returns:
                If pair.is_is(): the smallest stratum value S such that pair.distance <= S
                Otherwise, an empty string
            If strata not defined, returns None
        """
        if self.cis_strata and pair.is_cis:
            index = bisect.bisect_left(self.cis_strata, pair.distance)
            result = self.cis_strata[index]
            return result
        else:
            return None

    def classify(self, record: PairsSegment):
        """Generate a tuple with the classification code.
        
        The code can access the PairsSegment 'record', the stratum to which it
        maps 'stratum' (an integer if it is cis and strata are defined, an empty string
        if it is not cis and strata are defined, None if strata are undefined),
        and 'strata' (a list of integers or None).
        """
        
        if not self.compiled_classification_code:
            message = "No compiled classification code"
            if not self.conjuncts:
                message += "and list of conjuncts is empty."
            else:
                message += f"but list of conjuncts is {self.conjuncts}"
            raise TypeError(message)
        if not isinstance(record, PairsSegment):
            raise TypeError(f"Expected object of type PairsSegment, but got {type(record)}")
        
        
        strata = self.cis_strata
        stratum = self.get_stratum(record)
        result = eval(self.compiled_classification_code)
        return result
    
    def to_polars(self, distribution: DiscreteDistribution):
        """Output conjuncts plus 'count' as columns, rows as events + observed count
        
        Note: doesn't output the strata.
        """
        schema = self.conjuncts + ["count"]
        rows = [event + (count,) for event, count in distribution.items()]
        return DataFrame(rows, schema=schema, orient='row')
    
    def from_polars(self, df, cis_strata = None):
        """Treats all columns except 'count' as conjuncts, the 'count' column as the counts for each outcome
        """
        self.cis_strata = cis_strata
        self.conjuncts = [col for col in df.columns if col != "count"]
        self.__post_init__()
        distribution = DiscreteDistribution()
        for row in df.iter_rows():
            event = row[:-1]
            count = row[-1]
            distribution[event] = count
        return distribution

    def __getstate__(self):
        # Return a dictionary excluding the compiled code
        state = self.__dict__.copy()
        del state['compiled_classification_code']
        return state

    def __setstate__(self, state):
        # Restore the object's state and recompile the code
        self.__dict__.update(state)
        self.compiled_classification_code = compile(self.classification_code, '<string>', 'eval')

