import polars as pl

class DFView:
    def __init__(self, df):
        self.df = df

    def assert_conversion_allowed(self, fmt):
        conversions = type(self).conversions
        err = f"Conversion defined only for {conversions}, received {fmt}"
        assert fmt in conversions, err
    