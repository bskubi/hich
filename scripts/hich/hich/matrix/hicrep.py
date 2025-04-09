from hicrep.hicrep import sccByDiag
from hicrep.utils import *
import scipy.sparse as sp
import hictkpy
import numpy as np
import warnings
from pathlib import Path
from dataclasses import dataclass
import dataclasses
from typing import *

@dataclass
class HiCRepCaller:
    path1: Path | str
    path2: Path | str
    chromosome: str
    normalization: str
    resolution: int
    h: int
    d_bp_max: int
    b_downsample: int
    scc: float = None

    @property
    def files(self) -> Tuple[hictkpy.File, hictkpy.File]:
        """Return hictkpy.File handles for the two paths"""
        return (
            hictkpy.File(self.path1, self.resolution),
            hictkpy.File(self.path2, self.resolution)
        )

    @classmethod
    def call_scc(cls, caller: "HiCRepCaller"):
        file1, file2 = caller.files
        try:
            caller.scc = hicrep_scc_chrom(file1, file2, caller.chromosome, caller. normalization, caller.h, caller.d_bp_max, caller.b_downsample)
        except Exception as e:
            caller.scc = None
            warnings.warn(f"Exception on {caller}:\n{e}")

        return caller

    def as_dict(self):
        return {field.name: getattr(self, field.name) for field in dataclasses.fields(self)}


def hicrep_scc_chrom(file1: hictkpy.File, file2: hictkpy.File, chrom: str, normalization: str | None,
              h: int, dBPMax: int, bDownSample: bool) -> float:
    """Compute hicrep SCC score between two input hictkpy Files on a specific chromosome

    Args:
        file1: `hictkpy.File` Handle to the first file in the comparison
        file2: `hictkpy.File` Handle to the second file in the comparison
        chrom: `str`: The chromosome name to compare
        h: `int` Half-size of the mean filter used to smooth the
        input matrics
        dBPMax `int` Only include contacts that are at most this genomic
        distance (bp) away
        bDownSample: `bool` Down sample the input with more contacts
        to the same number of contacts as in the other input

    Returns:
        `float` scc score for the chromosome
    """
    # !Warning: this method has no specific unit test as of 2024/10/20 - Ben Skubi
    
    # Get resolution for each contact matrix
    resolution1 = file1.resolution()
    resolution2 = file2.resolution()

    # Input validation:
    # Same resolution
    # Both have the chromosome
    # Same number of bins per chromosome
    error_preamble = f"Input files {file1.path()} and {file2.path()}"

    assert resolution1 == resolution2,\
        f"{error_preamble} have different resolutions {resolution1} and {resolution2}"
    assert chrom in file1.chromosomes(), f"{file1.path()} does not contain chromosome {chrom}"
    assert chrom in file2.chromosomes(), f"{file2.path()} does not contain chromosome {chrom}"

    bins1 = file1.bins().query(f"chrom == '{chrom}'")
    bins2 = file2.bins().query(f"chrom == '{chrom}'")
    resolution = resolution1

    if resolution is None:
        # sometimes bin size can be None, e.g., input file has
        # non-uniform size bins.
        assert np.all(bins1[:] == bins2[:]),\
            f"Input contact matrix files don't have a unique bin size most likely "\
            f"because non-uniform bin size was used and the bins are defined "\
            f"differently for the two input contact matrix files"
        # In that case, use the median bin size
        resolution = int(np.median((bins1[:]["end"] - bins1[:]["start"]).values))
        warnings.warn(f"Input contact matrix files don't have a unique bin size most "\
                      f"likely because non-uniform bin size was used. HicRep "\
                      f"will use median bin size from the first contact matrix file "\
                      f"to determine maximal diagonal index to include", RuntimeWarning)

    assert len(bins1) == len(bins2),\
        f"{error_preamble} have different number of bins for chromosome {chrom} {len(bins1)} and {len(bins2)}"
    
    # Validate dBPMax
    if dBPMax == -1 or dBPMax is None:
        # this is the exclusive upper bound
        dMax = resolution
    else:
        dMax = dBPMax // resolution + 1

    assert dMax > 1, f"Input dBPmax is smaller than binSize"

    normalization = "NONE" if normalization is None else normalization

    # Get pixels
    p1 = file1.fetch(chrom).to_df()
    p2 = file2.fetch(chrom).to_df()
    file1.fetch()

    m1 = pixels2Coo(p1, bins1)
    m2 = pixels2Coo(p2, bins2)

    del p1
    del p2
    
    # get the total number of contacts as normalizing constant
    
    n1 = np.sum(m1.data[~np.isnan(m1.data)])
    n2 = np.sum(m2.data[~np.isnan(m2.data)])

    # Input validation on shape of chromosome matrix
    assert m1.size > 0, f"{error_preamble}, contact matrix 1 of chromosome {chrom} is empty"
    assert m1.shape[0] == m1.shape[1], f"{error_preamble}, contact matrix 1 of chromosome {chrom} is not square"

    assert m2.size > 0, f"{error_preamble}, contact matrix 2 of chromosome {chrom} is empty"
    assert m2.shape[0] == m2.shape[1], f"{error_preamble}, contact matrix 2 of chromosome {chrom} is not square"

    assert m1.shape == m2.shape, f"{error_preamble}, contact matrices of chromosome {chrom} have different input shape"

    return computeSCC(m1, m2, n1, n2, dMax, h, bDownSample)

def computeSCC(mS1: sp.coo_matrix,
               mS2: sp.coo_matrix,
               n1: float,
               n2: float,
               dMax: int,
               h: float,
               bDownSample: bool) -> float:

    """computeSCC 
    """
    # Compute scc score
    nDiags = mS1.shape[0] if dMax < 0 else min(dMax, mS1.shape[0])

    # remove major diagonal and all the diagonals >= nDiags
    # to save computation time
    m1 = trimDiags(mS1, nDiags, False)
    m2 = trimDiags(mS2, nDiags, False)
    del mS1
    del mS2

    if bDownSample:
        # do downsampling
        size1 = m1.sum()
        size2 = m2.sum()
        if size1 > size2:
            m1 = resample(m1, size2).astype(float)
        elif size2 > size1:
            m2 = resample(m2, size1).astype(float)
    else:
        # just normalize by total contacts
        m1 = m1.astype(float) / n1
        m2 = m2.astype(float) / n2

    if h > 0:
        # apply smoothing
        m1 = meanFilterSparse(m1, h)
        m2 = meanFilterSparse(m2, h)

    return sccByDiag(m1, m2, nDiags)
