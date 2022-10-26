#!/bin/bash
# set -eo pipefail

echo "Converting AIRR tab to nucleotide fasta ..."
python AIRRtab_to_fasta.py -i AIRR_sequences.tsv \
-o seqs.fasta

echo "Calling igblastn ..."
igblastn -germline_db_V databases/igblast_base/database/imgt_human_ig_v \
-germline_db_J databases/igblast_base/database/imgt_human_ig_j \
-germline_db_D databases/igblast_base/database/imgt_human_ig_d \
-organism human \
-query seqs.fasta \
-auxiliary_data databases/igblast_base/optional_file/human_gl.aux \
-show_translation \
-outfmt 19 > igblast_results.tab

echo "Converting igblastn result to protein fasta ..."
python igblast_to_fasta.py -i igblast_results.tab
