from multiprocessing import Pool
from typing import Tuple
from polars import DataFrame

def compute_pairs_stats_on_path(data):
    import time
    from hich.stats import DiscreteDistribution
    from hich.pairs import PairsClassifier, PairsFile, PairsSegment

    classifier, pairs_path = data
    pairs_file = PairsFile(pairs_path)
    stats = DiscreteDistribution()
    duration = 0

    for record in pairs_file:        
        outcome = classifier.classify(record)
        stats[outcome] += 1
    result = (pairs_path, stats)
    return result

def aggregate_classifier(pairs_stats_paths: list[str]) -> Tuple["PairsClassifier", list["DiscreteDistribution"]]:
    from hich.pairs import PairsClassifier
    from hich.stats import DiscreteDistribution
    dfs = []
    conjuncts = None
    raw_strata = set()
    for path in pairs_stats_paths:
        df = DataFrame(pairs_stats_header_tsv_path, separator = "\t", infer_schema_length = None)
        dfs.append(df)
        new_conjuncts = [col for col in df.columns if col != "count"]
        conjuncts = conjuncts or new_conjuncts
        assert conjuncts == new_conjuncts, "Conjuncts do not match for all hich stats files."
        if "stratum" in conjuncts:
            strata = set(df["stratum"].unique().to_list())
            raw_strata.update(strata)
    raw_strata = list(raw_strata)
    cis_strata = []
    for stratum in raw_strata:
        if stratum.isdigit():
            cis_strata.append(int(stratum))
        else:
            try:
                cis_strata.append(float(stratum))
            except ValueError:
                cis_strata.append("")
    classifier = PairsClassifier(conjuncts, cis_strata)
    distributions = []
    for df in dfs:
        distributions.append(classifier.from_polars(df))
    return classifier, distributions



def load_stats_and_classifier_from_file(pairs_stats_header_tsv_path):
    from hich.stats import DiscreteDistribution
    from hich.pairs import PairsClassifier, PairsFile
    df = DataFrame(pairs_stats_header_tsv_path, separator = "\t", infer_schema_length = None)
    conjuncts = [col for col in df.columns if col != "count"]
    cis_strata = None
    if "stratum" in conjuncts:
        raw_strata = df["stratum"].unique().to_list()
        cis_strata = []
        for s in raw_strata:
            if s.isdigit():
                cis_strata.append(int(s))
            elif s == "inf":
                cis_strata.append(float(s))
            else:
                cis_strata.append(s)
    classifier = PairsClassifier(conjuncts, cis_strata)
    distribution = DiscreteDistribution()
    for row in df.iter_rows(named = True):
        count = row["count"]
        outcome = (row[conjunct] for conjunct in conjuncts)
        distribution[outcome] = count
    return distribution, classifier

def compute_pairs_stats_on_path_list(classifier, pairs_paths):
    data = [(classifier, path) for path in pairs_paths]
    return Pool().map(compute_pairs_stats_on_path, data)