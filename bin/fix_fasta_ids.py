import argparse
import os.path
from Bio import SeqIO

def parse():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--fasta_f",
        help = "file containing processed results"
    )
    args = parser.parse_args()
    fix_fasta(**vars(args))

def fix_fasta(fasta_f):
    records = SeqIO.parse(fasta_f, 'fasta')
    outf = ((os.path.basename(fasta_f)).split("."))[0] + ".adj.fasta"
    rec_list = []
    processed = []
    for record in records:
        desc = record.description
        annot_split = desc.split(";")
        annot_sub = ([x for x in annot_split if 'locus_tag' in x])[0]
        locus_tag = annot_sub.rsplit('locus_tag=', 1)[1]
        if locus_tag not in processed:
            processed.append(locus_tag)
            record.id = locus_tag
            record.description = ''
            rec_list.append(record)
    with open(outf, "w") as output_handle:
        SeqIO.write(rec_list, output_handle, "fasta")

if __name__ == "__main__":
    parse()
