from dataclasses import dataclass, field
from smart_open import smart_open
from hich.parse.pairs_header import PairsHeader
from hich.parse.pairs_segment import PairsSegment
from pathlib import Path

@dataclass
class PairsFile:
    filepath_or_object: str = None
    mode: str = None
    header: PairsHeader = field(default_factory = PairsHeader)

    def __init__(self,
        filepath_or_object,
        mode = "rt",
        template = None,
        text = None,
        header = None,
        check_header = True,
        check_segments = True):
        self.mode = mode
        if filepath_or_object:
            file_did_not_exist = not Path(filepath_or_object).exists()
            self.filepath_or_object = smart_open(filepath_or_object, mode = self.mode)
            if "r" in self.mode:
                self.header = PairsHeader()
                lines = []
                while (line := self.filepath_or_object.readline()).startswith("#"):
                    lines.append(line)
                    if line.startswith("#columns:"):
                        break
                header_text = "".join(lines)
                self.header = PairsHeader.from_text(header_text)
            elif "w" in self.mode or ("a" in self.mode and file_did_not_exist):
                self.header = header
                self.filepath_or_object.write(self.header.to_string())

    def __del__(self):
        self.close()

    def close(self):
        self.filepath_or_object.close()

    def pair_segment_from_text(self, line):
        stripped = line.strip()
        if not stripped:
            raise StopIteration
        fields = stripped.split()
        field_vals = {dict(enumerate(self.header.columns))[idx]: val for idx, val in enumerate(fields)}
        return PairsSegment(**field_vals)

    def __iter__(self):
        return self

    def __next__(self):
        line = self.filepath_or_object.readline()
        return self.pair_segment_from_text(line)

    def to_header(self):
        self.filepath_or_object.seek(0)
    
    def to_records(self, record_number=0):
        self.to_header()

        # First, scan for the first non-comment line
        while True:
            position = self.filepath_or_object.tell()  # Get the current file pointer position
            line = self.filepath_or_object.readline()  # Read a line manually

            if not line:  # If we reach EOF, exit the loop
                break

            if not line.startswith("#"):
                self.filepath_or_object.seek(position)  # Seek back to the start of the non-comment line
                if record_number == 0:
                    return line

        # If record_number is not 0, continue reading lines
        while record_number > 0:
            line = self.filepath_or_object.readline()
            record_number -= 1
            if not line:  # Handle the case where EOF is reached
                return None
            if record_number == 0:
                return line


    def write(self, pairs_segment: PairsSegment):
        self.filepath_or_object.write("\n" + pairs_segment.to_string())