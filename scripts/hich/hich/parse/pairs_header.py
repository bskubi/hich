from pathlib import Path

class PairsHeader:
    version_prefix = "## pairs format v"

    def __init__(self, chromsizes = {}, columns = [], command = []):
        self.chromsizes = chromsizes
        self.columns = columns
        self.command = command

    def to_dict(self):
        return self.__dict__
    
    def to_string(self):
        version_line = PairsHeader.version_prefix + str(self.pairs_format_version)
        chromsizes_lines = [
            f"#chromsize: {chrom} {size}"
            for chrom, size
            in self.chromsizes.items()
        ]
        columns_line = "#columns: " + " ".join(self.columns)
        other_lines = []
        for key, val in self.to_dict().items():

            if key not in ["pairs_format_version", "chromsizes", "columns"]:
                if isinstance(val, list):
                    for elem in val:
                        other_lines.append(f"#{key}: {elem}")
                else:
                    other_lines.append(f"#{key}: {val}")
        header_lines = [version_line] + chromsizes_lines + other_lines + [columns_line]
        return "\n".join(header_lines)

    @classmethod
    def from_dict(self, header_dict):
        header = PairsHeader()
        header.update(header_dict)
    
    @classmethod
    def from_text(self, from_text):
        header = PairsHeader()
        lines = from_text.split("\n")
        assert lines[0].startswith(PairsHeader.version_prefix), f"Pairs must start with ## pairs format v1.0 but first line was {line}"
        header.pairs_format_version = lines[0].removeprefix("## pairs format v")
        for line in lines[1:]:
            if not line.startswith("#"):
                break
            fields = line.split()
            field_type = fields[0]
            rest = line.removeprefix(field_type).lstrip()
            if field_type == "#chromsize:":
                contig, size = fields[1:]
                header.chromsizes[contig] = size
            elif field_type == "#command:":
                header.command.append(rest)
            elif field_type == "#columns:":
                header.columns = fields[1:]
            elif field_type.endswith(":"):
                field_name = field_type[1:-1]
                if field_name not in self.__dict__:
                    setattr(header, field_name, rest)
                elif isinstance(header.__dict__[field_name], str):
                    header.__dict__[field_name] = [header.__dict__[field_name], rest]
                else:
                    header.__dict__[field_name].append(rest)
        return header
    
    def __repr__(self):
        return repr(str(self))

    def __str__(self):
        return self.to_string()