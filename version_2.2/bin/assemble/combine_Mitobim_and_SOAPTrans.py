#!/usr/bin/python3
"""
combine_Mitobim_and_SOAPTrans.py

Copyright (c) 2017-2018 Guanliang Meng <mengguanliang@foxmail.com>.

This file is part of MitoZ.

MitoZ is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MitoZ is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MitoZ.  If not, see <http://www.gnu.org/licenses/>.

"""

#import sys
import argparse

description = """
combine the sequences of Mitobim and SOAPTrans, i.e. connect the sequence of
SOAPTrans and the terminal regions (extension part) of Mitobim. Only if the
mapping happens will the Mitobim extension part be merged with SOAPTrans results.
"""

parser = argparse.ArgumentParser(description=description)
parser.add_argument("-trans", required=True, metavar="<STR>",
					help="SOAPTrans fasta file")
parser.add_argument("-mitobim", required=True, metavar="<STR>",
					help="Mitobim fasta file")
parser.add_argument("-len", type=int, default=100, metavar="[INT]",
						help="alignment length default: %(default)s")
parser.add_argument("-mis", type=int, default=5, metavar="[INT]",
					help="mismatch allowed for alignment default: %(default)s")
parser.add_argument("-o", required=True, metavar="<STR>",
			help="outfile (seperated sequences from SOAPTrans and Mitobim)")
parser.add_argument("-o2", required=True, metavar="<STR>",
					help="outfile (Only the combined sequence)")

args = parser.parse_args()

def readfa2dict(infile):
	seq_dict = {}
	seq_id = ""
	seq = ""
	first = 1
	fh_in = open(infile, 'r')
	for line in fh_in:
		line = line.rstrip()
		if line.startswith(">"):
			if first == 1:
				first = 0
				seq_id = line.split()[0].replace("_bb", "")
				continue
			#seq = seq.upper()
			seq_dict[seq_id] = seq
			seq_id = line.split()[0].replace("_bb", "")
			seq = ""
		else:
			seq += line
	seq_dict[seq_id] = seq
	fh_in.close()
	return seq_dict

def find_index(str1, anotherseq, mismatch_allowed):
	# length of str1 and str2 must be equal
	length = len(str1)
	j = 0
	for j in range(0, len(anotherseq)-length+1):
		str2 = anotherseq[j:j+length]
		i = 0
		mismatch = 0
		for i in range(0, len(str1)):
			if str1[i] != str2[i]:
				mismatch += 1
		if mismatch <= mismatch_allowed:
			#print("dfjakdfjkj")
			return j
	return -1

def vertical_bar_by_identity(seq1, seq2):
	# seq1 and seq2 must have equal length
	pass


trans_dict = readfa2dict(args.trans)
mitobim_dict = readfa2dict(args.mitobim)


fh_out = open(args.o, 'w')
fh_out2 = open(args.o2, 'w')
for seq_id in trans_dict.keys():
	trans_seq = trans_dict[seq_id]
	mitobim_seq = mitobim_dict[seq_id]

	trans_subseq_5 = trans_seq[0:args.len].lower()
	#print("trans 5' 10bp:", trans_subseq_5)
	trans_subseq_3 = trans_seq[-args.len:].lower()
	#print("trans 3' 10bp:", trans_subseq_3)

	message = "and this end of Mitobim sequence will not be merged with " + \
			"SOAPTrans sequence"

	index_5 = find_index(trans_subseq_5, mitobim_seq, args.mis)
	index_3 = find_index(trans_subseq_3, mitobim_seq, args.mis)

	combined_seq = ""
	if (index_5 != -1) and (index_3 != -1):
		combined_seq = mitobim_seq[0:index_5] + trans_seq + mitobim_seq[(index_3+args.len):]
	elif (index_5 == -1) and (index_3 != -1):
		combined_seq = trans_seq + mitobim_seq[(index_3+args.len):]
	elif (index_5 != -1) and (index_3 == -1):
		combined_seq = mitobim_seq[0:index_5] + trans_seq
	else:
		print(seq_id, trans_seq, sep="\n", file=fh_out)
		print(seq_id, mitobim_seq, sep="\n", file=fh_out)
		print(seq_id, trans_seq, sep="\n", file=fh_out2)
		continue

	if combined_seq:
		print(seq_id, trans_seq, sep="\n", file=fh_out)
		print(seq_id, combined_seq, sep="\n", file=fh_out)
		print(seq_id, mitobim_seq, sep="\n", file=fh_out)
		print(seq_id, combined_seq, sep="\n", file=fh_out2)

fh_out.close()




