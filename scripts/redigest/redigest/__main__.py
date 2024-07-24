import click
from Bio import SeqIO
from Bio.Restriction import RestrictionBatch
from Bio.Seq import Seq
import sys
from itertools import pairwise
from smart_open import smart_open

def get_cut_sites(seq_record, enzymes):
    cut_sites = []
    analysis = enzymes.search(seq_record.seq)
    for enzyme in analysis:
        positions = [0] + analysis[enzyme]
        for pos in pairwise(positions):
            cut_sites.append((seq_record.id, *pos))
    if len(cut_sites) > 0:
        cut_sites.append((seq_record.id, cut_sites[-1][2], len(seq_record.seq)))
    else:
        cut_sites = [(seq_record.id, 0, len(seq_record.seq))]
    return cut_sites

def write_bed_file(cut_sites, output_file):
    with open(output_file, "w") as bed_file:
        for site in cut_sites:
            bed_file.write(f"{site[0]}\t{site[1]}\t{site[2]}\n")

def main(reference_genome, enzyme_names, output_file):
    # Load the reference genome
    handle = smart_open(reference_genome, "rt")
    seq_record = SeqIO.parse(handle, "fasta")
    genome = SeqIO.to_dict(seq_record)

    enzymes = RestrictionBatch(enzyme_names)
    
    # Process each sequence record in the genome
    cut_sites = []
    for seq_record in genome.values():
        cuts = get_cut_sites(seq_record, enzymes)
        cut_sites.extend(cuts)
    
    # Write to BED file
    write_bed_file(cut_sites, output_file)

@click.command()
@click.option("--output", default = None, help = "Output file")
@click.argument("reference")
@click.argument("enzymes", nargs = -1)
def redigest(output, reference, enzymes):
    enzymes = set(enzymes)
    if "Arima" in enzymes:
        enzymes.add("DpnII")
        enzymes.add("HinfI")
        enzymes.remove("Arima")
    main(reference, enzymes, output)

if __name__ == "__main__":
    redigest()