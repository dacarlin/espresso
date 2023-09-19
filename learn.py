import json 
from itertools import product 
import argparse 
from collections import Counter 

from biotite.sequence.io.fasta import FastaFile 

from espresso.lib import CodonSequence 


parser = argparse.ArgumentParser()
parser.add_argument("fasta", help="FASTA file containing coding sequences") 
args = parser.parse_args() 
fasta = FastaFile.read(args.fasta)

codons = []
for header, record in fasta.items():
    codon_sequence = CodonSequence(record) 
    for codon in codon_sequence.codons:
        codons.append(codon) 


result = dict(Counter(codons))

print(json.dumps(result, indent=4))
