#!/usr/bin/env python3

import argparse
import magnumopus

parser = argparse.ArgumentParser(description="Perform in-silico PCR on two assemblies and align the amplicons")
parser.add_argument("-1", dest="ASSEMBLY1", type=str, required=True, help="Path to the first assembly file")
parser.add_argument("-2", dest="ASSEMBLY2", type=str, required=True, help="Path to the second assembly file")
parser.add_argument("-p", dest="PRIMERS", type=str, required=True, help="Path to the primer file")
parser.add_argument("-m", dest="MAX_AMPLICON_SIZE", type=int, required=True, help="maximum amplicon size for isPCR")
parser.add_argument("--match", dest="MATCH", type=int, required=True, help="match score to use in alignment")
parser.add_argument("--mismatch", dest="MISMATCH", type=int, required=True, help="mismatch penalty to use in alignment")
parser.add_argument("--gap", dest="GAP", type=int, required=True, help="gap penalty to use in alignment")
args = parser.parse_args()
ASSEMBLY1 = args.ASSEMBLY1
ASSEMBLY2 = args.ASSEMBLY2
PRIMERS = args.PRIMERS
MAX_AMPLICON_SIZE = args.MAX_AMPLICON_SIZE
MATCH = args.MATCH
MISMATCH = args.MISMATCH
GAP = args.GAP

pcr1 = magnumopus.ispcr(PRIMERS, ASSEMBLY1, MAX_AMPLICON_SIZE)
seq1 = pcr1.split('\n')[1]
pcr2 = magnumopus.ispcr(PRIMERS, ASSEMBLY2, MAX_AMPLICON_SIZE)
seq2 = pcr2.split('\n')[1]
complement = {"A" : "T", "C": "G", "G" : "C", "T" : "A"}
seq2_comp = "".join(complement[base] for base in seq2)
seq2_revcomp = seq2_comp[::-1]
(aln1, aln2), score = magnumopus.needleman_wunsch(seq1, seq2, MATCH, MISMATCH, GAP)
(aln1_revcomp, aln2_revcomp), score_revcomp = magnumopus.needleman_wunsch(seq1, seq2_revcomp, MATCH, MISMATCH, GAP)

if score > score_revcomp:
	print(aln1)
	print(aln2)
	print(score)
else:
	print(aln1_revcomp)
	print(aln2_revcomp)
	print(score_revcomp)