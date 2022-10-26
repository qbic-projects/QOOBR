import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", type=str, help="Path to input AIRR format tab-separated (TSV) table coming from igblast.")
parser.add_argument("-o", "--outname", type=str, default="sequences_aa.fasta", help="Name of output fasta file.")

args = parser.parse_args()

tab = pd.read_csv(args.input, sep="\t")
tab["fasta_desc"] = tab["v_call"] + " " + tab["d_call"] + " " + tab["j_call"]

out_file = args.outname
ids = tab["sequence_id"]
descriptions = tab["fasta_desc"]
seqs_nt = tab["sequence"]
seqs_aligned_aa = tab["sequence_alignment_aa"]

with open(args.outname, 'wt') as out_handle:
    for (ii, desc, seq_nt, seq_al_aa) in zip(ids, descriptions, seqs_nt, seqs_aligned_aa):

        seq_orig_nt = Seq(seq_nt)

        # Translating sequence
        seq_translated = seq_orig_nt.translate()

        # If there is a stop codon, move to next frame
        if '*' in seq_translated:
            seq_translated = seq_orig_nt[1:].translate()
            if '*' in seq_translated:
                seq_translated = seq_orig_nt[2:].translate()
                # If all frames have stop codons, fall back to aligned aa sequence from igblast
                if '*' in seq_translated:
                    print("Falling back to igblast for sequence {}".format(ii))
                    seq_translated = Seq(seq_al_aa)

        record = SeqRecord(
            seq_translated,
            id = ii,
            name = ii,
            description = str(desc),
        )
        SeqIO.write(record, out_handle, "fasta")
