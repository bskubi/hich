from hich.parse.pairs_file import PairsFile
from hich.coverage.trans_cis_thresholds import TransCisThresholds
from hich.coverage.pairs_histogram import PairsHistogram

pairs_file = PairsFile("240802_Arima_reseq_dedup.pairs.gz")
cis_thresholds = [10, 20, 50, 100]
space = TransCisThresholds(ignore_code = "not pair.ur()", cis_thresholds = cis_thresholds)

dist = PairsHistogram(space)

for i, pair in enumerate(pairs_file):
    dist.count_pair(pair)
    if i >= 1000:
        break

print(dist)
print(dist.to_probdist())
prob = [.7, .15, .05, .1, 0]
print(dist.downsample_to_probdist(prob))
print(dist.downsample_to_probdist(prob).to_count(50))

"""
0. On input
    Associate samples with downsampling groups
        Default: all in their own target profile group

    Specify a minimum number of reads per contig.
        - Above: follow normal methods
        - Below: either retain all or drop all

1. To form target probdist PairHistograms, convert count PairHistograms to
probabilities and average them, leaving profiles unchanged if they are alone
in the sample group.
    - Can also provide an explicit target for each sample group
    - Specify some sample files and/or contigs to ignore

2. The user can specify:
    2a) (optional) Downsample to target probdist for sample group
    2b) (optional) Downsample to any of:
        - Specific size
        - Specific fraction
        - Minimum sample size within sample group
        - Minimum sample size across sample groups

3. This yields new PairsHistograms for each sample. We then use MultiSampler
to selection sample the pairs files.
"""