from typing import Optional, Annotated, NewType
from pydantic import AfterValidator

def valid_resolutions(bin_resolutions: list[int]) -> list[int]:
    if any([r <= 0 for r in bin_resolutions]):
        raise ValueError(f"Contact matrix binning resolutions {bin_resolutions} must be all positive integers.")
    min_r = min(bin_resolutions)
    if any([r % min_r != 0 for r in bin_resolutions]):
        raise ValueError(f"Contact matrix binning resolutions {bin_resolutions} must all be multiples of minimum value: {min_r}.")
    return bin_resolutions

MatrixBinResolutionsValidator = AfterValidator(valid_resolutions)

DEFAULT_BIN_RESOLUTIONS = [
    10, 20, 50,
    100, 200, 500,
    1_000, 2_000, 5_000, 
    10_000, 20_000, 50_000, 
    100_000, 200_000, 500_000, 
    1_000_000, 2_000_000, 5_000_000
]