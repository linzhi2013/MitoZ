#!/usr/bin/python3
"""
gc_and_coverage_of_fasta.py

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

from Bio import SeqIO
from Bio.Graphics.GenomeDiagram import Diagram
from reportlab.lib.colors import yellow, red, green, orange
import argparse

from gc_content import gc_content

parser = argparse.ArgumentParser(description="")
parser.add_argument("-i", metavar="<STR>", required=True, help="input fasta file")
parser.add_argument("-d", metavar="<STR>", required=True, help="file of sequencing depth along the sequences. format per-line: seq_id x1 x2 ... xn")
parser.add_argument("-o", metavar="<STR>", required=True, help="output prefix")
parser.add_argument("-w", metavar="<INT>", type=int, default=100, help="window size for caculating GC content along the sequences. default: %(default)s")

args = parser.parse_args()

fasta_file = args.i
seq_depth_file = args.d
prefix = args.o
window_size = args.w

## coverage or depth
cov_dict = {}
fh_in = open(seq_depth_file, 'r')
for i in fh_in:
        i = i.strip()
        line = i.split()
        seqid = line[0]
        depth = [(j-1,int(line[j])) for j in range(1,len(line)-1)]
        cov_dict[seqid] = depth
fh_in.close()


for record in SeqIO.parse(fasta_file, 'fasta'):
	name = record.id
	outfile = prefix + ".%s" % name


	gdd_circular = Diagram("mitogenome")
	#gdd_linear = Diagram("mitogenome")

	## GC content
	gdt_circular_gc = gdd_circular.new_track(2, name="GC Content", greytrack=True, greytrack_labels=1, start=0, end=len(record), scale=True, scale_fontangle=45, scale_format="SInt", scale_largetick_interval=1000, scale_smallticks=0.05, scale_smalltick_interval=100, scale_ticks=True, axis_labels=True)
	gdgs_circular_gc = gdt_circular_gc.new_set('graph')
	gc_list = gc_content(str(record.seq), window_size)[0]
	gdgs_circular_gc.new_graph(gc_list, "GC content", style="line")

	## Coverage
	gdt_circular_depth = gdd_circular.new_track(3, name="Coverage", greytrack=True, greytrack_labels=1, start=0, end=len(record), scale=True, scale_format="SInt", scale_largetick_interval=1000, scale_smallticks=0.05, scale_smalltick_interval=100, scale_ticks=True, axis_labels=True)
	gdgs_circular_depth = gdt_circular_depth.new_set('graph')
	if record.id in cov_dict:
		depth_list = cov_dict[record.id]
		gdgs_circular_depth.new_graph(depth_list, "Coverage", style="bar") # color=colors.lightgreen, altcolor=colors.red
	else:
		print("can not find coverage data of %s" % record.id)
		sys.exit(0)

	gdd_circular.draw(format="circle", circular=False, orientation="landscape", pagesize='A4', fragments=1, start=0, end=len(record))
	gdd_circular.write(outfile + ".coverage-GC_circle.pdf", "PDF")
	gdd_circular.write(outfile + ".coverage-GC_circle.svg", "SVG")






