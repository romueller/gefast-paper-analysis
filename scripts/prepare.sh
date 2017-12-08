#!/bin/bash

# Code for (un)even data set adapted from:
# Mahé F, Rognes T, Quince C, de Vargas C, Dunthorn M. (2014) 
# Swarm: robust and fast clustering method for amplicon-based studies. 
# PeerJ 2:e593 https://doi.org/10.7717/peerj.593
#
# Supplement 8: https://doi.org/10.7717/peerj.593/supp-8
#
# Downloads and preprocesses the data for the subsequent analyses



SWARM=tools/Swarm-2.1.13

# ===== ELDERMET data set =====

# 1) Download
# see http://www.ebi.ac.uk/ena/data/view/SRP003158
runs134=(134392 134393 134394 134396)
runs135=(135535 135536 135537 135538)
runs136=(136479 136480 136481 136482 136483 136484 136485 136486 136487 136488 136489 136490 136491)

for r in ${runs134[@]}
do
	wget -P data/eldermet ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR134/SRR${r}/SRR${r}.fastq.gz
done
for r in ${runs135[@]}
do
	wget -P data/eldermet ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR135/SRR${r}/SRR${r}.fastq.gz
done
for r in ${runs136[@]}
do
	wget -P data/eldermet ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR136/SRR${r}/SRR${r}.fastq.gz
done

# 2) Combine into a single FASTA file
cat data/eldermet/SRR*.fastq.gz > data/eldermet/eldermet.fastq.gz
seqtk seq -a data/eldermet/eldermet.fastq.gz > data/eldermet/eldermet_raw.fasta
rm data/eldermet/*.fastq.gz

# 3) Remove sequences containing n/N
python scripts/removeN.py data/eldermet/eldermet_raw.fasta data/eldermet/eldermet.fasta

# 4) Dereplicate
./${SWARM} -d 0 -a 1 -w data/eldermet/eldermet_derep.fasta data/eldermet/eldermet.fasta > /dev/null


# ===== Even data set =====

# 1) Download
wget -P data/even http://sbr2.sb-roscoff.fr/download/externe/de/fmahe/even.fasta.bz2
bzip2 -d data/even/even.fasta.bz2
mv data/even/even.fasta data/even/even_raw.fasta

# 2) Remove sequences containing n/N
python scripts/removeN.py data/even/even_raw.fasta data/even/even.fasta

# 3) Dereplicate
./${SWARM} -d 0 -a 1 -w data/even/even_derep.fasta data/even/even.fasta > /dev/null


# ===== Uneven data set =====

# 1) Download
wget -P data/uneven http://sbr2.sb-roscoff.fr/download/externe/de/fmahe/uneven.fasta.bz2
bzip2 -d data/uneven/uneven.fasta.bz2
mv data/uneven/uneven.fasta data/uneven/uneven_raw.fasta

# 2) Remove sequences containing n/N
python scripts/removeN.py data/uneven/uneven_raw.fasta data/uneven/uneven.fasta

# 3) Dereplicate
./${SWARM} -d 0 -a 1 -w data/uneven/uneven_derep.fasta data/uneven/uneven.fasta > /dev/null


