import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", type=str, help="Path to input AIRR format tab-separated (TSV) table.")
parser.add_argument("-o", "--outname", type=str, default="sequences.fasta", help="Name of output fasta file.")

args = parser.parse_args()

tab = pd.read_csv(args.input, sep="\t", low_memory=False)
tab['unique_sequence_id'] = tab.apply(lambda x: f'{x["sequence_id"]}_{x["sample_id"]}', axis=1)

# check that unique_sequence_id is really unique
assert len(tab) == len(tab['unique_sequence_id'].unique())

out_file = args.outname
ids = tab["unique_sequence_id"]
seqs = tab["sequence"]
with open(args.outname, 'wt') as out_handle:
    for (ii, seq) in zip(ids, seqs):
        record = SeqRecord(
            Seq(seq),
            id=ii,
            name=ii,
            description='',
        )
        SeqIO.write(record, out_handle, "fasta")