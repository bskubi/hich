from hich.pairs import PairsFile, PairsHeader, PairsSegment
from hypothesis import given, example, strategies as st, assume, target, settings
from tests.pairs import *
from io import StringIO

# Data can be downloaded from https://data.4dnucleome.org/files-processed/4DNFIRMZ7QTE/
# Should be saved to tests/pairs/data
