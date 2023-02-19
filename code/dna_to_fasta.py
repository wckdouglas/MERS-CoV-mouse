from Bio import SeqIO
from pathlib import Path
import click

def parse_snapgene(dna_file):
    with dna_file.open('rb') as dna:
        records = SeqIO.parse(dna, "snapgene")

        for record in records:
            print(f">{record.id} {dna_file.name}\n{record.seq.upper()}")

@click.command()
@click.option('--data-path', required=True, help='data directory to search for .dna files')
def main(data_path: str):
    """
    Searching snapgene file in the data_path and convert to fasta

    :param str data_path: data directory to search for .dna files
    """
    for dna_file in Path(data_path).glob("*.dna"):
        parse_snapgene(dna_file)

if __name__ == "__main__":
    main()

