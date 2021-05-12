#!/usr/bin/python3
"""
mitogenome_visualizaion_from_genbank.py

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

import sys
from Bio import SeqIO
#from Bio.Graphics.GenomeDiagram import Diagram
from GenomeDiagram import Diagram
from reportlab.lib.colors import yellow, red, green, orange
import argparse
#from gc_content_along_the_sequence import GC_along_the_sequence
#from gc_content import gc_content

description = """
draw a linear and a circular SVG file for input mitochondrial genome in Genbank file.
"""
parser = argparse.ArgumentParser(description=description)
parser.add_argument("-i", metavar="<STR>", required=True, help="genbank infile")
parser.add_argument("-p", metavar="<STR>", default="", help="prefix of outfiles when not use '-o'. default: species name of input genbank file")
parser.add_argument("-o", metavar="<STR>", help="outfile prefix. For drawing all diagrams in this one outfile when mutiple sequences in genbank infile.")
parser.add_argument("-c", action="store_true", default=False, help="all sequences are cicular? default: %(default)s")

args = parser.parse_args()

genbankfile = args.i
circular = args.c
#prefix = args.p

if args.o:
	gdd = Diagram("mitogenome")

	#tracks_dict = {}

	max_len = 0
	j = 0
	for record in SeqIO.parse(genbankfile, 'gb'):

		seq_len = len(record)
		max_len = max(max_len, seq_len)

		#tracks_dict.setdefault(record.id, [])
		#tracks_dict[record.id].append(seq_len)
		#tracks_dict[record.id].append(track_num)

		track_num = 4 + j

		#name = record.annotations['organism'].replace(" Unclassified.", "")
		#outfile = name.replace(" ", "_")
		#outfile = name + "_" + record.id
		#outfile = outfile.replace(" ", "_")

		gdt = gdd.new_track(track_num, name=record.id, greytrack=True, greytrack_labels=1, start=0, end=len(record), scale=True, scale_format="SInt", scale_largetick_interval=1000, scale_smallticks=0.1, scale_smalltick_interval=100, scale_ticks=True)
		gdfs = gdt.new_set()
		j += 3

		for fea in record.features:
			if fea.type == "CDS":
				gdfs.add_feature(fea, sigil="ARROW", color=yellow, label=True, label_position="middle")
			elif fea.type == "tRNA":
				gdfs.add_feature(fea, sigil="ARROW", color=green, label=True, label_position="middle")
			elif fea.type == "rRNA":
				gdfs.add_feature(fea, sigil="ARROW", color=red, label=True, label_position="middle")

	if circular:
		gdd.draw(format="circle", circular=circular,  pagesize='A4', fragments=1, start=0, end=max_len)
	else:
		gdd.draw(format="circle", circular=circular,  pagesize='A4', fragments=1, start=0, end=max_len+20)
	gdd.write(args.o + "_circle.pdf", "PDF")
	gdd.write(args.o + "_circle.svg", "SVG")

	gdd.draw(format="linear", circular=circular,  pagesize='A2', fragments=1, start=0, end=max_len)
	gdd.write(args.o + "_linear.pdf", "PDF")
	gdd.write(args.o + "_linear.svg", "SVG")

else:
	for record in SeqIO.parse(genbankfile, 'gb'):
		gdd_circular = Diagram("mitogenome")
		gdd_linear = Diagram("mitogenome")

		name = record.annotations['organism'].replace(" Unclassified.", "")
		#outfile = name.replace(" ", "_")
		if args.p != "":
			outfile = args.p + "_" + record.id
		else:
			outfile = name + "_" + record.id
		outfile = outfile.replace(" ", "_")

		## new_track(), scale 是否画中间的黑线
		gdt_circular = gdd_circular.new_track(3, name=name, greytrack=True, greytrack_labels=0, start=0, end=len(record), scale=True, scale_format="SInt", scale_largetick_interval=1000, scale_smallticks=0.1, scale_smalltick_interval=100, scale_ticks=True)
		gdt_linear = gdd_linear.new_track(3, name=name, greytrack=True, greytrack_labels=0, start=0, end=len(record), scale=True, scale_format="SInt", scale_largetick_interval=1000, scale_smallticks=0.1, scale_smalltick_interval=100, scale_ticks=True)

		gdfs_cicular = gdt_circular.new_set()
		gdfs_linear = gdt_linear.new_set()


		## circular
		for fea in record.features:
			if fea.type == "CDS":
				gdfs_cicular.add_feature(fea, sigil="ARROW", color=yellow, label=True, label_position="middle")
			elif fea.type == "tRNA":
				gdfs_cicular.add_feature(fea, sigil="ARROW", color=green, label=True, label_position="middle")
			elif fea.type == "rRNA":
				gdfs_cicular.add_feature(fea, sigil="ARROW", color=red, label=True, label_position="middle")
			elif fea.type in ["D-loop", "misc_feature", "rep_origin", "repeat_region"]:
				gdfs_cicular.add_feature(fea, sigil="ARROW", color=orange, label=True, label_position="middle")

		## linear
		for fea in record.features:
			if fea.location.strand == 1:
				if fea.type == "CDS":
					gdfs_linear.add_feature(fea, sigil="ARROW", color=yellow, label=True, label_position="middle")
				elif fea.type == "tRNA":
					gdfs_linear.add_feature(fea, sigil="ARROW", color=green, label=True, label_position="middle")
				elif fea.type == "rRNA":
					gdfs_linear.add_feature(fea, sigil="ARROW", color=red, label=True, label_position="middle")
				elif fea.type in ["D-loop", "misc_feature", "rep_origin", "repeat_region"]:
					gdfs_linear.add_feature(fea, sigil="ARROW", color=orange, label=True, label_position="middle")
			else:
				if fea.type == "CDS":
					gdfs_linear.add_feature(fea, sigil="ARROW", color=yellow, label=True, label_angle=120, label_position="middle")
				elif fea.type == "tRNA":
					gdfs_linear.add_feature(fea, sigil="ARROW", color=green, label=True, label_angle=120, label_position="middle")
				elif fea.type == "rRNA":
					gdfs_linear.add_feature(fea, sigil="ARROW", color=red, label=True, label_angle=120, label_position="middle")
				elif fea.type in ["D-loop", "misc_feature", "rep_origin", "repeat_region"]:
					gdfs_linear.add_feature(fea, sigil="ARROW", color=orange, label=True, label_angle=120, label_position="middle")

		if circular:
			gdd_circular.draw(format="circle", circular=circular,
				orientation="landscape", pagesize='A4', fragments=1,
				start=0, end=len(record))
		else:
			gdd_circular.draw(format="circle", circular=circular,
				orientation="landscape", pagesize='A4',
				fragments=1, start=0, end=len(record)+20)
		gdd_circular.write(outfile + ".gene_circle.pdf", "PDF")
		gdd_circular.write(outfile + ".gene_circle.svg", "SVG")

		gdd_linear.draw(format="linear", circular=circular, orientation="landscape", pagesize="A2", fragments=1, start=0, end=len(record))
		gdd_linear.write(outfile + ".gene_linear.pdf", "PDF")
		gdd_linear.write(outfile + ".gene_linear.svg", "SVG")







