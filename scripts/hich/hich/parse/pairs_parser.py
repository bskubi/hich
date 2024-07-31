import polars as pl
import smart_open
import warnings
from hich.parse.pbgzip import *


class PairsParser:
    def __init__(self, filename):
        self.filename = filename
        self.write_file = None
        
    def columns_row(self):
        """Returns the 0-indexed offset of the final header row with the column names

        row_index: 0-indexed offset of the final header row (header rows start with '#')
        columns_dict: {colname:index} dict of which columns correspond to which column names

        If a row starting with '#columns:' is present it is assumed to be the final header row.
        If no '#columns:' row is present, returns the last header row with the minimal seven
        column names given by the .pairs specification.
        
        Returns: (row_index, columns_dict)
        """

        default_columns = ['readID', 'chrom1', 'pos1', 'chrom2', 'pos2', 'strand1', 'strand2']
        row = 0

        with smart_open.open(self.filename, "rt", encoding='utf-8') as file:
            for line in file:
                if line.startswith('#columns:'):
                    column_names = line.split()[1:]
                    return (row, column_names)
                elif line.startswith('#'):
                    row += 1
                else:
                    return (row-1, default_columns)
        return (row - 1, default_columns)

    def header(self):
        """Returns list of all lines starting with '#' until but not including the first line not starting with '#'
        """
        header = []
        with smart_open.open(self.filename, "rt", encoding="utf-8") as file:
            for line in file:
                if line.startswith('#'):
                    header.append(line)
                else:
                    break
        return header

    def batch_iter(self, n_rows):
        def read_chunk(skip_rows, columns):
            return pl.read_csv(self.filename,
                             skip_rows = skip_rows,
                             n_rows = n_rows,
                             raise_if_empty = False,
                             separator = '\t',
                             has_header = False,
                             new_columns = columns,
                             dtypes = {"readID":pl.String,
                                       "chrom1":pl.String,
                                       "chrom2":pl.String})
            
        columns_row, columns = self.columns_row()

        skip_rows = columns_row + 1
        with open(self.filename) as file:
            df = read_chunk(skip_rows, columns)
            skip_rows += len(df)
            yield df
            while len(df) == n_rows:
                df = read_chunk(skip_rows, columns)
                skip_rows += len(df)
                if not df.is_empty():
                    yield df

    def write_append(self, filename, df = None, header_end = None):
        warnings.filterwarnings("ignore", message="Polars found a filename")

        if self.write_file is None:
            
            self.write_file = smart_open.open(filename, "w", compression = compression(filename))
            header_no_columns = self.header()[:-1]
            final_header_line_fields = header_no_columns[-1].split()
            
            if header_end:
                pn_field = [field[3:]
                    for field in final_header_line_fields
                    if field.startswith("PN:")]
                
                header_end.PP = pn_field[0] if pn_field else "null"
                
            header_no_columns = "".join(header_no_columns)
            columns = " ".join(["#columns:"] + df.columns) + "\n"
            header_end = str(header_end)
            header = header_no_columns + header_end + columns
            self.write_file.write(header)

        if df is not None:
            df.write_csv(self.write_file,
                        include_header = False,
                        separator = '\t')

    def close(self):
        if self.write_file is not None:
            filename = self.write_file.name
            self.write_file.close()
            self.write_file = None
            if filename.endswith(".gz") or filename.endswith(".gzip"):
                pbgzip_compress(filename)