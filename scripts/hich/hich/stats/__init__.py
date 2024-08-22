from collections import Counter
from hich.parse.pairs_file import PairsFile
from numpy import searchsorted

def pair_stats(pairs_file: PairsFile, output = None, conjuncts = list[str], cis_strata = None, print = False):
    events = Counter()

    cols = ', '.join(conjuncts)
    event_code = f"({cols})"
    event_code_compiled = compile(event_code, '<string>', 'eval')

    nat_partition = lambda cuts: sorted(list(set(cuts + [float('inf')]))) if cuts else None
    strata = nat_partition(cis_strata)

    if strata:
        get_stratum = lambda pair: strata[searchsorted(strata, pair.distance)] if pair.is_cis() else ""
    else:
        get_stratum = lambda pair: None

    for i, pair in enumerate(pairs_file):
        stratum = get_stratum(pair)
        event = eval(event_code_compiled)
        events[event] += 1

    return events