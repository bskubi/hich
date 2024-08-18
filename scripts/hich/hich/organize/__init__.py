import polars as pl
import pysam
import Bio
from hich.parse.seqio_splitter import SeqIOSplitter
from hich.parse.alignment_splitter import AlignmentSplitter
from hich.parse.pairs_file import PairsFile
from hich.parse.pairs_segment import PairsSegment
from hich.parse.pairs_splitter import PairsSplitter
from pathlib import Path

def annotation_dict(filename, has_header, separator):
    df = pl.read_csv(filename, separator = separator, has_header = has_header)
    col1, col2 = df.columns[:2]
    return dict(zip(df[col1], df[col2]))

def organize_seqio_paired(fmt: str,
                          f1: str,
                          f2: str,
                          out_dir: str,
                          annotations: str,
                          head: int,
                          record_code: str,
                          filename_code: str):
    handles = zip(Bio.SeqIO.parse(f1), Bio.SeqIO.parse(f2))
    splitters = [SeqIOSplitter(fmt), SeqIOSplitter(fmt)]
    for i, ends in enumerate(handles):
        if i >= head:
            break

        record1, record2 = ends
        filename1, filename2 = eval(filename_code)
        exec(record_code)
        splitters[0].write(filename1, record1)
        splitters[1].write(filename2, record2)

def organize_seqio_single(fmt: str,
                          f1: str,
                          out_dir: str,
                          annotations: str,
                          head: int,
                          record_code: str,
                          filename_code: str):
    handle = Bio.SeqIO.parse(f1)
    splitter = SeqIOSplitter(fmt)
    for i, record in enumerate(handle):
        if i >= head:
            break

        filename = eval(filename_code)
        exec(record_code)
        splitter.write(filename, record)


def organize_alignment(filename: str,
                       out_dir: str,
                       annotations: str,
                       head: int,
                       record_code: str,
                       filename_code: str):

    with pysam.AlignmentFile(filename) as alignments:
        splitter = AlignmentSplitter(header = alignments.header)
        for i, record in enumerate(alignments):
            if head and i == head:
                break

            filename = eval(filename_code, {'record':record, 'annotations':annotations})
            if record_code: exec(record_code, {'record': record, 'annotations': annotations})
            path = Path(out_dir) / Path(filename)
            
            splitter.write(str(path), record)

def organize_pairs(f1: str,
                   annotations: str,
                   head: str,
                   record_code: str,
                   filename_code: str):
    with PairsFile(filename) as pairs:
        splitter = PairsSplitter(header = pairs.header)
        for i, record in enumerate(pairs):
            if head and i == head:
                break

            filename = eval(filename_code, {'record':record, 'annotations':annotations})
            if record_code: exec(record_code, {'record': record, 'annotations': annotations})
            path = Path(out_dir) / Path(filename)
            
            splitter.write(str(path), record)

def organize(fmt: str,
             f1: str,
             f2: str = None,
             out_dir: str = None,
             annot_file: str = None,
             annot_has_header: str = False,
             annot_separator: str = "\t",
             head: int = None,
             record_code: str = None,
             filename_code: str = None):
    annotations = annotation_dict(annot_file, annot_has_header, annot_separator)
    out_dir = out_dir or ""

    if record_code: record_code = compile(record_code, '<string>', 'exec')
    if filename_code: filename_code = compile(filename_code, '<string>', 'eval')
    
    if fmt in ["fastq", "fasta"]:
        if f2: organize_seqio_paired(fmt, f1, f2, out_dir, annotations, head, record_code, filename_code)
        else: organize_seqio_single(fmt, f1, out_dir, annotations, head, record_code, filename_code)
    elif fmt in ["sam", "bam"]:
        organize_alignment(f1, out_dir, annotations, head, record_code, filename_code)
    elif fmt in ["pairs"]:
        organize_pairs(f1, out_dir, annotations, head, record_code, filename_code)