"""Parse 4DN .pairs format headers
https://github.com/4dn-dcic/pairix/blob/master/pairs_format_specification.md

Works with compressed files using smart_open. The pairtools suite saves using pbgzip, which smart_open can't use by default. The smart_open_with_pbgzip makes both smart_open and DuckDB use pbgzip.
"""

import smart_open_with_pbgzip
from smart_open import smart_open
from pathlib import Path
from dataclasses import dataclass
from typing import *

@dataclass
class PairsHeader:
    """Parser for 4DN .pairs format header

    https://github.com/4dn-dcic/pairix/blob/master/pairs_format_specification.md

    path: path to the .pairs file
    lines: list of raw header lines
    columns: list of column names'
    misformatted: true if most recently parsed header was misformatted
    data: extracted from non-columns header lines
    
    Methods:
        __init__ -- calls parse by default
        parse -- extract header from file or raw text

    Properties and setters
        parse_header_to_columns -- convert #columns: line in header to self.columns
        parse_columns_to_header -- convert self.columns to #columns: line in header
        text -- raw text of header
        non_columns_text -- raw text of header, any #columns: line removed
        non_columns_lines -- raw text of header in list of separate lines (has setter, does not reparse)
        last_line_columns -- returns True if last header line starts with #columns:, false otherwise 
    """
    path: Path | None = None
    lines: List[str] | None = None
    columns: List[str] = None
    misformatted: bool = None
    data: Dict[str, Any] = dict

    def __init__(self, 
                 data: Path | str, 
                 from_text: bool = False, 
                 columns: List[str] | None = None, 
                 parse: bool = True, 
                 parse_path: bool | None = None, 
                 replace_columns: bool | None = True,
                 replace_data: bool = True
        ):
        """
        Arguments:
            data: Either raw header text (parsed this way if from_text is true) or a path to a 4DN .pairs file (from_text is false)
            from_text: Determines whether to parse data as raw text or path
            columns: List of predetermined columns. Either appended to lines (final line is not a columns line) or replaces final columns line
            parse: Immediately parse the file or text 
            parse_path: If parse is called, enforces whether it parses self.path or self.lines. If None, uses what's available, raises exception if both/None are
            replace_columns: If columns is set, then if replace_columns is True, columns is replaced if a #columns: line is in the header
        """
        if from_text:
            # rstrip the line because a trailing \n on the last line will be set by split as a trailing empty line 
            self.lines = [line + "\n" for line in data.rstrip().split("\n")]
            
            if columns:
                columns_line = " ".join(["#columns:", *columns]) + "\n"
                if self.lines[-1].startswith("#columns:"):
                    self.lines[-1] = columns_line
                else:
                    self.lines.append(columns_line)
        else:
            self.path = Path(data)
        
        if (self.lines or self.path.exists()) and parse:
            self.parse(parse_path, replace_columns, replace_data)

    def parse(self, parse_path: bool = None, replace_columns: bool = False, replace_data: bool = False):
        """Parse 4DN .pairs header or raw text lines

        If parsing .pairs file, file may be compressed or non-compressed. Works with compressed files using smart_open. The pairtools suite saves using pbgzip, which smart_open can't use by default. The smart_open_with_pbgzip makes both smart_open and DuckDB use pbgzip.
        
        If misformatted, PairsHeader.misformatted is set to True.
        Extracts columns to PairsHeader.columns.
        Stores other lines as strings, lists or a dict (for "chromsize") in data.

        Arguments:
            store_lines: Store the raw lines of the header in PairsHeader.header
            parse_path: If true, parses path. If false, parse lines. If None, parse whichever exists (raise if both exist).
            replace_columns: If columns is set, then if replace_columns is True, columns is replaced if a #columns: line is in the header
            replace_data: If data set set, then if replace_data is True, old data is cleared before update.

        File format description:
        https://github.com/4dn-dcic/pairix/blob/master/pairs_format_specification.md
        """
        parse_path = self.update_parse_path_HELPER(parse_path)
        update_columns = (self.columns and replace_columns) or not self.columns
        
        error_preamble = self.get_error_preamble_HELPER(parse_path)
        header = self.get_header_HELPER(parse_path, error_preamble)
        self.columns_HELPER(header, update_columns)
        self.data_HELPER(header, replace_data)

    ###################################
    #   Getters and setters
    ###################################

    @property
    def text(self):
        """Get text of all lines"""
        return "".join(self.lines) if self.lines else None
    
    @property
    def non_columns_text(self):
        """Get text of all lines except columns, if present"""
        return "".join(self.non_columns_lines) if self.lines else None
    
    @property
    def last_line_columns(self) -> bool:
        """Returns whether the last header line starts with #columns:"""
        return self.lines[-1].startswith("#columns:") if self.lines else False

    @property
    def columns_text(self):
        """Get text of columns line"""
        if self.columns:
            joined_columns = " ".join(["#columns:", *self.columns])
        if self.lines[-1].startswith("#columns:"):
            columns_text = self.lines[-1]
            assert not self.columns or columns_text == joined_columns, f"Both self.columns and columns line were set, but did not match.\nColumns line: {columns_line}.\nself.columns: {self.columns}"
        elif self.columns:
            columns_text = joined_columns
        else:
            columns_text = None
        return columns_text

    @property
    def non_columns_lines(self):
        """Get lines except final #columns: line"""
        non_columns_lines = None
        if self.lines:
            non_columns_lines = self.lines[:-1] if self.last_line_columns else self.lines
        return non_columns_lines

    @non_columns_lines.setter
    def non_columns_lines(self, non_columns_lines: List[str]):
        """Set lines other than #columns:"""
        if self.last_line_columns:
            self.lines = non_columns_lines + self.lines[-1]
        else:
            self.lines = non_columns_lines

    def parse_header_to_columns(self):
        """Set self.columns based on #columns: line in header"""
        assert self.last_line_columns, f"Tried to parse header to columns but last line does not start with '#columns:'\nLast line: {self.lines[-1]}"
        self.columns = self.lines[-1].split()[1:]

    def parse_columns_to_header(self):
        """Append/replace #columns: line in header based on self.columns"""
        assert self.columns, f"Tried to parse columns to header but self.columns is empty or None"
        new_columns_line = " ".join(["#columns:", *self.columns])
        if self.last_line_columns:
            self.columns[-1] = new_columns_line
        else:
            self.columns.append(new_columns_line)

    #####################################
    # Helper methods
    #####################################

    def update_parse_path_HELPER(self, parse_path: bool):
        """Set parse_path, if None, to true/false based on whether self.path or self.lines is available"""
        if parse_path is None:
            if self.path and not self.lines:
                parse_path = True
            elif self.lines and not self.path:
                parse_path = False
            elif self.path and self.lines:
                raise ValueError("Called parse on PairsHeader with parse_path set to None, which is ambiguous. Set to True or False when both path and lines are set before call to parse.")
        if not (self.path or self.lines):
            raise ValueError("Called parse on PairsHeader where neither path nor lines are available, so there is nothing to parse.")
        return parse_path
    
    def load_pairs_header_HELPER(self, path: Path | str) -> List[str]:
        """Load list of lines from 4DN .pairs file."""
        path = Path(path)
        error_preamble = f"Attempted to parse pairs header from {path} with absolute path {path.absolute()}"
        assert path.exists(), f"{error_preamble}, but file does not exist."
        self.lines = []
        try:
            with smart_open(path, "rt", encoding="utf-8") as f:
                for line in f:
                    if line.startswith("#"):
                        self.lines.append(line)
                    else:
                        break
        except Exception as e:
            print(f"{error_preamble}, but caught an exception while opening and parsing the file: {e}")
            raise

    def get_header_HELPER(self, parse_path: bool, error_preamble: str):
        """Load header from file or self.lines"""
        if parse_path:
            self.load_pairs_header_HELPER(self.path)
        header = self.lines
        if not header:
            self.misformatted = True
            raise ValueError(f"{error_preamble}, but data did not begin with a header where lines start with #.")
        return header

    def get_error_preamble_HELPER(self, parse_path: bool):
        """Get start of error message depending on if parsing from path or text lines"""
        if parse_path:
            error_preamble = f"Attempted to parse pairs header from {self.path} with absolute path {self.path.absolute()}"
        else:
            error_preamble = f"Attempted to parse pairs header from lines"
        return error_preamble

    def columns_HELPER(self, header: List[str], update_columns: bool):
        """Build self.columns from header"""
        if self.last_line_columns and update_columns:
            self.parse_header_to_columns()
        if not self.columns:
            self.misformatted = True

    def data_HELPER(self, header: List[str], replace_data: bool):
        """Build self.data dict from header"""
        if replace_data:
            self.data = {}

        for line in header[:-1]:
            trimmed = line.lstrip("#").strip()
            if not trimmed:
                continue
            split_line = line.split()

            # Trim leading '#' and trailing ":"
            field = split_line[0]
            value = split_line[1:]

            if field == "chromsize":
                self.data.setdefault("chromsize", {})
                chrom, size = value
                self.data["chromsize"][chrom] = size
            else:
                entry = " ".join(value)
                if field in self.data:
                    if isinstance(self.data[field], str):
                        self.data[field] = [self.data[field]]
                    self.data[field].append(entry)
                else:
                    self.data[field] = entry