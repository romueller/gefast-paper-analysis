# Adapted from:
# https://www.biostars.org/p/183279/
#
# Removes sequences containing n/N



import sys
from Bio import SeqIO

handle = open(sys.argv[1], "rU")
output_handle = open(sys.argv[2], "w")

filtered = [record for record in SeqIO.parse(handle, "fasta") if ((record.seq.count('N') == 0) and (record.seq.count('n') == 0))]
SeqIO.write(filtered, output_handle, "fasta")

output_handle.close()
handle.close()
