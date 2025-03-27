import click
import hich.fasta.digest_re
import polars
import smart_open_with_pbgzip
from smart_open import smart_open

@click.group()
def fasta():
    pass

@fasta.group()
def digest():
    pass

@digest.command
@click.argument("in_path")
@click.argument("out_path")
@click.argument("enzymes_and_protocols", nargs=-1)
def re(in_path, out_path, enzymes_and_protocols):
    """Create BED file of restriction fragments

    \b
    IN_PATH: FASTA file to digest
    OUT_PATH: Name of BED file to output digest breakpoints
    ENZYMES_AND_PROTOCOLS: One or more space-separated REBASE restriction enzyme names (case sensitive), or protocols or Hi-C kit names (case-insensitive).

    \b
    Examples:
        hich fasta digest-re hg38.fasta hg38_digest.bed DpnII DdeI
        hich fasta digest-re hg38.fasta hg38_digest.bed hic3.0
        hich fasta digest-re hg38.fasta hg38_digest.bed Arima

    \b        
    In addition to any combination of REBASE enzymes, you can invoke hich fasta digest re with any of the following protocol and kit names:

    \b
    ENZYMES_AND_PROTOCOLS arguments and corresponding enzymes:
    "HiC 3.0" --> Enzymes: ['DpnII', 'DdeI']
    "Hi-C 3.0" --> Enzymes: ['DpnII', 'DdeI']
    "Arima Genome-Wide HiC+" --> Enzymes: ['DpnII', 'HinfI']
    "Arima" --> Enzymes: ['DpnII', 'HinfI']
    "Phase Proximo 2021+ Plant" --> Enzymes: ['DpnII']
    "Phase Plant" --> Enzymes: ['DpnII']
    "Phase Proximo 2021+ Animal" --> Enzymes: ['DpnII']
    "Phase Animal" --> Enzymes: ['DpnII']
    "Phase Proximo 2021+ Microbiome" --> Enzymes: ['Sau3AI', 'MluCI']
    "Phase Microbiome" --> Enzymes: ['Sau3AI', 'MluCI']
    "Phase Proximo 2021+ Human" --> Enzymes: ['DpnII']
    "Phase Human" --> Enzymes: ['DpnII']
    "Phase Proximo 2021+ Fungal" --> Enzymes: ['DpnII']
    "Phase Fungal" --> Enzymes: ['DpnII']

    \b
    Notes:
    1. Protocol names are case-insensitive, so you can call hich fasta digest_re [IN_PATH] [OUT_PATH] 'hic 3.0'.
    2. Use quotes for protocol names with spaces.
    3. The enzyme list is filtered for uniqueness so repeated enzyme names don't affect results.
    """
    handle = smart_open(out_path, "wt")
    for cutsites in hich.fasta.digest_re.digest_re(in_path, enzymes_and_protocols):
        for row in cutsites.iter_rows():
            outputs = [str(it) for it in row]
            text = "\t".join(outputs) + "\n"
            handle.write(text)
