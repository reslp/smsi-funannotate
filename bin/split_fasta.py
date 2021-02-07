#!/usr/bin/env python
# written by Philipp Resl
import os
import sys
from Bio import SeqIO
import argparse

if sys.version_info[0] < 3:
    raise Exception("Must be using Python 3")

pars = argparse.ArgumentParser(prog="split_fasta.py", description = """This script will split a fasta file into multiple files.""", epilog = """written by Philipp Resl""")
pars.add_argument('--fasta', dest="fasta_file", required=True, help="Path to FASTA file")
pars.add_argument('--n_records', dest="n_records", required=True, help="Number of records in each file")
pars.add_argument('--outdir', dest="outdir", required=True, help="Path to output directory.")
args=pars.parse_args()

#make paths absolute
if not os.path.exists(args.outdir):
	os.makedirs(args.outdir)
if not os.path.isabs(args.outdir):
	args.outdir = os.path.abspath(args.outdir)

if not os.path.isabs(args.fasta_file):
	args.fasta = os.path.abspath(args.fasta_file)

def write_sequences(chunk, nchunk):
	filename = os.path.basename(args.fasta_file)
	filename += "_"+str(nchunk)
	file = open(args.outdir+"/"+filename, "w")
	file.write(chunk)
	file.close()

i = 0
chunk = ""
nchunk = 0
for record in SeqIO.parse(args.fasta, "fasta"):
	if i < int(args.n_records):
		chunk += ">"+str(record.id) + "\n"
		chunk += str(record.seq) + "\n"
		i += 1
	else:
		write_sequences(chunk, nchunk)
		i = 0
		chunk = ""	
		nchunk += 1
