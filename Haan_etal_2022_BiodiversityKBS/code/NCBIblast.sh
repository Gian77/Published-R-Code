#!/bin/bash 

# Copyright © 2022 Gian M.N. Benucci, Ph.D.
# email: benucci@msu.edu URL: https://github.com/Gian77

# To run this script you need to download the whole NCBI nucleotide locally
# Usage: sh NCBIblast.sh.sh sequence.fasta <organism_group>

export BLASTDB="/mnt/home/benucci/DATABASES/NCBI_1121/"

file="$1"
organism="$2"

conda activate blast_tools

blastn \
    -task megablast \
    -num_threads $SLURM_CPUS_PER_TASK \
    -db /mnt/home/benucci/DATABASES/NCBI_1121/nt \
    -query $file \
    -max_target_seqs 10 \
    -max_hsps 10 \
    -outfmt '6 qseqid staxids bitscore pident evalue length qlen slen qcovs qcovhsp sskingdoms scomnames sscinames sblastnames stitle' \
    -evalue 1e-5 \
    -perc_identity 60.00 \
    -out results_${organism}.txt 

# extract BLAST results id-taxid-bitscore-pident-evalue-coverage
grep "OTU" results_${organism}.txt | sed '/OTU/s/\t/\|/g' | cut -f 1,2,3,4,5,9 -d "|" | sed '/OTU/s/|/\t/g' > res_${organism}.txt

# extract taxids and generate taxonomy lineage form taxids
cut -f 1,2 results_${organism}.txt | cut -f 1 -d ";" > taxids_${organism}.txt

grep "OTU" taxids_${organism}.txt | cut -f 2 | taxonkit lineage -t --data-dir /mnt/home/benucci/DATABASES/NCBI_1121/ | sed 's/\t/\;/g' | cut -f 1,3,5,6,7,8,9,10 -d ";" | sed 's/;/\t/g' > taxonomy_${organism}.txt

paste res_${organism}.txt taxids_${organism}.txt taxonomy_${organism}.txt > taxformat_${organism}.txt

sed '1i seqID\ttaxID\tbitscore\tp-identity\te-value\tq-coverage\tseqID\ttaxID\ttaxID\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies' taxformat_${organism}.txt  > table_taxformat_${organism}.txt 

conda deactivate
