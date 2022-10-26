import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", type=str, help="Path to input AIRR format tab-separated (TSV) table.")
parser.add_argument("-o", "--outname", type=str, default="sequences.fasta", help="Name of output fasta file.")

args = parser.parse_args()

tab = pd.read_csv(args.input, sep="\t")
tab["fasta_desc"] = tab["v_call_genotyped"] + " " + tab["d_call"] + " " + tab["j_call"]

out_file = args.outname
ids = tab["sequence_id"]
descriptions = tab["fasta_desc"]
seqs = tab["sequence"]
with open(args.outname, 'wt') as out_handle:
    for (ii, desc, seq) in zip(ids, descriptions, seqs):
        record = SeqRecord(
            Seq(seq),
            id=ii,
            name=ii,
            description=str(desc),
        )
        SeqIO.write(record, out_handle, "fasta")