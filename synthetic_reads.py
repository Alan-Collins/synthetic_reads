#!/usr/bin/env python

import sys
import os
import argparse
import random
import pathlib
import numpy as np



class ReadPair():
	def __init__(self, id, seq1='', seq2='', qual1='', qual2=''):
		self.rname = f'read_{id}'
		self.seq = {1: seq1, 2: seq2}
		self.qual = {1: qual1, 2: qual2}

	@staticmethod
	def rand_qual(read_len, mean, sd):
		mean += 33 # adjust phred
		qual = []
		for i in range(read_len):
			qual_score = int(np.random.normal(mean, sd))
			if qual_score > 75:
				qual_score = 75
			if qual_score < 33:
				qual_score = 33
			qual.append(chr(qual_score))

		return ''.join(qual)


	@staticmethod
	def rand_seq(n):
		seq = []
		for i in range(n):
			seq.append(random.choice(['A', 'T', 'C', 'G']))

		return ''.join(seq)


	@staticmethod
	def rev_comp(string):
		"""Reverse complement a string of nucleotide sequence
		
		Args:
		string (str):
			Nucleotide sequence

		Returns:
		str:
			reverse complement of given nucleotide sequence
		
		Raises:
		TypeError: If string is not str.

		"""
		if type(string) is not str:
			raise TypeError(
				"string must be str, not {}.".format(type(string).__name__))

		rev_str = ''
		rev_comp_lookup = {
		"A" : "T", 
		"T" : "A", 
		"C" : "G", 
		"G" : "C", 
		"a" : "t", 
		"t" : "a", 
		"c" : "g", 
		"g" : "c",
		}
		for i in reversed(string):
			if i in "ATCGatcg":
				rev_str += rev_comp_lookup[i]
			else:
				rev_str += i
		return rev_str

	@staticmethod
	def sequence_fragment(seq, read_len, fragment_size):
		mid = random.randint(0,len(seq))
		fragment_len = random.randint(fragment_size[0], fragment_size[1])
		start = int(mid-fragment_len/2)
		end = int(mid+fragment_len/2)
		
		if start < -(read_len-50):
			r1 = ReadPair.rand_seq(read_len)
		else:
			r1 = seq[max([start, 0]): start+read_len]
			if len(r1) < read_len:
				r1 = ReadPair.rand_seq(read_len-len(r1)) + r1
		
		if end > len(seq)+read_len-50:
			r2 = ReadPair.rand_seq(read_len)
		else:
			r2 = ReadPair.rev_comp(seq[end-read_len: min([len(seq), end])])
			if len(r2) < read_len:
				r2 += ReadPair.rand_seq(read_len-len(r2))

		switch = random.choice([True, False])
		if switch:
			r1, r2 = r2, r1
		
		return r1, r2


	def to_fastq_strings(self):
		read1 = (
			f'@{self.rname}\n'
			f'{self.seq[1]}\n'
			'+\n'
			f'{self.qual[1]}\n'
		)
		
		read2 = (
			f'@{self.rname}\n'
			f'{self.seq[2]}\n'
			'+\n'
			f'{self.qual[2]}\n'
		)

		return read1, read2

	@classmethod
	def from_seq(cls, seq, id, qual_mean, qual_sd, read_len, fragment_size):
		seq1, seq2 = ReadPair.sequence_fragment(seq, read_len, fragment_size)
		qual1 = ReadPair.rand_qual(read_len, qual_mean, qual_sd)
		qual2 = ReadPair.rand_qual(read_len, qual_mean, qual_sd)

		return cls(id, seq1, seq2, qual1, qual2)



class Reads():
	def __init__(self, reads_dict):
		self.reads = reads_dict

	def write_fastq(self, outprefix):
		reads1, reads2 = [], []
		for read_pair in self.reads.values():
			a, b = read_pair.to_fastq_strings()
			reads1.append(a)
			reads2.append(b)
		
		with open(outprefix + '1.fastq', 'w') as fout:
			fout.write("".join(reads1))
		
		with open(outprefix + '2.fastq', 'w') as fout:
			fout.write("".join(reads2))

		


	@classmethod
	def from_fasta(cls, fasta_dict, num_reads, qual_mean, qual_sd, read_len, fragment_size):
		reads = {}
		input_seqs = [s for s in fasta_dict.values()]
		for n in range(num_reads):
			seq = random.choice(input_seqs)
			reads[n] = ReadPair.from_seq(seq, n, qual_mean, qual_sd, read_len, fragment_size)

		return cls(reads)



def cmdline_args():

	p = argparse.ArgumentParser(
		description=""
		)
	p.add_argument(
		"-n", "--num_reads", required = False, type=int, default=500,
		help="number of reads to be generated (default: %(default)s)"
		)
	p.add_argument(
		"-q", "--qual_mean", required = False, type=int, default=30,
		help="average quality of bases in reads (default: %(default)s)"
		)
	p.add_argument(
		"-s", "--qual_sd", required = False, type=int, default=10,
		help="standard deviation of distribution of base qualities (default: %(default)s)"
		)
	p.add_argument(
		"-l", "--read_len", required = False, type=int, default=250,
		help="Length of reads (default: %(default)s)"
		)
	p.add_argument(
		"-i", "--insert_size_range", required = False, type=str, default="400:700",
		help="range of desired insert sizes. Format 'min:max' (default: %(default)s)"
		)
	p.add_argument(
		"-f", "--fasta", required = True, type=str,
		help="path to fasta sequence for which you want to generate reads"
		)
	p.add_argument(
		"-o", "--outprefix", required = False, type=str, default="./synthreads_",
		help="prefix/path to prepend to output files (default: %(default)s)"
		)
	p.add_argument(
		"--seed", required = False, type=int, default=None,
		help="Set the seed to control randomness of read generation"
		)

	return p.parse_args()


def fasta_to_dict(FASTA_file):
	"""Read a fasta file into a dict 

	Dict has headers (minus the > symbol) as keys and the associated 
	sequence as values.
	
	Args:
	  FASTA_file (str): 
		path to fasta format file

	Returns:
	  dict: 
		dict of format {fasta_header : sequence}

	Raises:
	  TypeError: If FASTA_file is not a str
	  OSError: If FASTA_file is not the path to an existing file
	"""
	
	if type(FASTA_file) is not str:
		raise TypeError(
			"FASTA_file must be str, not {}.".format(type(FASTA_file).__name__))

	if not os.path.exists(FASTA_file):
		raise OSError(
			"FASTA_file must be the path to an existing file.")


	fasta_dict = {}
	with open(FASTA_file, 'r') as f:
		multifasta = f.read()
	f.close()
	fastas = multifasta.split(">")
	trimmed_fastas = []
	for i in fastas:
		if len(i) != 0:
			trimmed_fastas.append(i)

	fastas = trimmed_fastas

	for i in fastas:
		header = i.split("\n")[0]
		seq = "".join(i.split("\n")[1:])
		fasta_dict[header] = seq

	return fasta_dict


def make_outdir(outprefix):
	# make outdir if doesn't exist
	if outprefix[-1] == "/":
		if not os.path.exists(outprefix):
			pathlib.Path(outprefix).mkdir(parents=True, exist_ok=True)
	else:
		directory = "/".join(outprefix.split("/")[:-1])
		if not os.path.exists(directory):
			pathlib.Path(directory).mkdir(parents=True, exist_ok=True)


def main(args):
	make_outdir(args.outprefix)
	random.seed(args.seed)
	fasta_dict = fasta_to_dict(args.fasta)
	reads = Reads.from_fasta(
		fasta_dict,
		args.num_reads,
		args.qual_mean,
		args.qual_sd,
		args.read_len,
		[int(i) for i in args.insert_size_range.split(":")]
	)
	reads.write_fastq(args.outprefix)

	
if __name__ == '__main__':
	args = cmdline_args()
	main(args)	
