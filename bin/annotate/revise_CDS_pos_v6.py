#!/usr/bin/python
"""
revise_CDS_pos_v6.py

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

usage = '''
Description
    Revising CDS position by translation.

Usage
    python3 revise_CDS_pos.py  <in.fasta>  <candidate.CDS.position3>  <table_num>  <outfile>

Version
    1.1
    update: when there are more than one same occurs in sam sequences, will output all results.

Author
    Guanliang MENG, BGI
    email: mengguanliang@genomics.cn
'''

## subroutine


## main
###############################################################################################

if len(sys.argv) != 5:
    print(usage)
    sys.exit(0)

in_f1, in_f2, table, out_f = sys.argv[1:5]

table = int(table)

TAA_note = "TAA stop codon is completed by the addition of 3' Aresidues to the mRNA"
pos_dict = {}
pos_revised_dict = {} ## [start_pos, end_pos, dire, TAA_note]
ref_dict = {}
stop_still_missing_dict = {}
start_still_missing_dict = {}

occur_times = {}

################################################################################################

fh = open(in_f2, 'r')
for line in fh:
    line = line.rstrip()
    seq_id, gene_name, ref_len, ref_start, ref_end, scaf_start, scaf_end, dire = line.split()
    ref_len = int(ref_len)
    ref_start = int(ref_start)
    ref_end = int(ref_end)
    scaf_start = int(scaf_start)
    scaf_end = int(scaf_end)

    if (seq_id in occur_times) and (gene_name in occur_times[seq_id]):
        # v6
        occur_times[seq_id][gene_name] += 1
        occur_time = occur_times[seq_id][gene_name]
        #occur_time = occur_times[seq_id][gene_name] + 1
        gene_name = gene_name + "-" + str(occur_time)
    else:
        if seq_id not in occur_times:
            occur_times.setdefault(seq_id, {})
        occur_times[seq_id][gene_name] = 1

    pos_revised_dict.setdefault(seq_id, {})[gene_name] = [-1,-1,dire,'NA']
    pos_dict.setdefault(seq_id, {})[gene_name] = (scaf_start, scaf_end, dire)
    ref_dict.setdefault(seq_id, {})[gene_name] = (ref_start, ref_end, ref_len)

fh.close()

###########################################

for record in SeqIO.parse(in_f1, 'fasta'):
    seq_id = record.id
    if seq_id not in pos_dict:
        continue
    seq_len = len(record)
    for gene_name in pos_dict[seq_id]:
        scaf_start, scaf_end, dire = pos_dict[seq_id][gene_name]
        #print('>',seq_id, gene_name, scaf_start, scaf_end, dire)

        #####################
        ## find stop codon ##
        start_pos = -1
        end_pos = -1
        sub_seq_stop = ''
        stop_codon_missing = 1
        j = 0
        while j < 10:
            if dire == '+':
                if (scaf_end - 1) + j*3 + 1 < len(record.seq):
                    start_pos = (scaf_end - 1) + 1
                    end_pos = (scaf_end - 1) + 1 + j*3
                    sub_seq_stop = record.seq[start_pos:end_pos]
                else:
                    break
            elif dire == '-':
                if (scaf_start-1) - j*3 >= 0:
                    start_pos = (scaf_start-1) - j*3
                    end_pos = (scaf_start - 1)
                    sub_seq_stop = record.seq[start_pos:end_pos]
                    sub_seq_stop = sub_seq_stop.reverse_complement()
                else:
                    break
            else:
                print("Error: %s direction is unknown!\n" % gene_name)
                sys.exit(0)

            pro_seq_stop = sub_seq_stop.translate(table=table)
            if pro_seq_stop.endswith('*'):
                if dire == '+':
                    pos_revised_dict[seq_id][gene_name][1] = end_pos
                    #print(end_pos)
                elif dire == '-':
                    #pos_revised_dict[seq_id][gene_name][1] = start_pos+1
                    pos_revised_dict[seq_id][gene_name][0] = start_pos+1
                    #print(start_pos+1)
                stop_codon_missing = 0
                break

            j += 1

        ## if stop codon not found, try to find 'T', 'TA'
        if stop_codon_missing:
            j = 0
            while j < 10:
                if dire == '+':
                    if (scaf_end - 1) + 1 + j*3 < len(record):
                        start_pos = (scaf_end - 1) + 1 + j*3
                        ## in case the DNA sequence is in lower case
                        sub_T = record.seq[start_pos:(start_pos+1)].upper()
                        sub_TA = record.seq[start_pos:(start_pos+2)].upper()
                        # must try 'TA' first, then try 'T'
                        if sub_TA == 'TA':
                            pos_revised_dict[seq_id][gene_name][1] = start_pos+1+1
                            pos_revised_dict[seq_id][gene_name][3] = TAA_note
                            stop_codon_missing = 0
                            break
                        elif sub_T == 'T':
                            pos_revised_dict[seq_id][gene_name][1] = start_pos+1
                            pos_revised_dict[seq_id][gene_name][3] = TAA_note
                            stop_codon_missing = 0
                            break
                    else:
                        break
                elif dire == '-':
                    if (scaf_start - 1) - 1 - j*3 >= 0:
                        start_pos = (scaf_start - 1) - 1 - j*3
                        sub_T = record.seq[start_pos:(start_pos+1)].upper()
                        sub_TA = record.seq[(start_pos-1):(start_pos+1)].upper()
                        sub_T = sub_T.reverse_complement()
                        sub_TA = sub_TA.reverse_complement()
                        # must try 'TA' first, then try 'T'
                        if sub_TA == 'TA':
                            pos_revised_dict[seq_id][gene_name][0] = start_pos-1+1
                            pos_revised_dict[seq_id][gene_name][3] = TAA_note
                            stop_codon_missing = 0
                            break
                        elif sub_T == 'T':
                            #print(start_pos+1, "TAA stop codon is completed by the addition of 3' Aresidues to the mRNA")
                            #pos_revised_dict[seq_id][gene_name][1] = start_pos+1
                            pos_revised_dict[seq_id][gene_name][0] = start_pos+1
                            pos_revised_dict[seq_id][gene_name][3] = TAA_note
                            stop_codon_missing = 0
                            break
                    else:
                        break

                j += 1

        ## if stop codon still missing, keep this gene_name and the original positions
        if stop_codon_missing:
            stop_still_missing_dict.setdefault(seq_id, []).append(gene_name)
            if dire == '+':
                pos_revised_dict[seq_id][gene_name][1] = ">" + str(scaf_end)
            elif dire == '-':
                #if start_pos == 0:
                #    start_pos = 2
                pos_revised_dict[seq_id][gene_name][0] = "<" + str(scaf_start)

        ## find stop codon ##
        #####################


        ######################
        ## find start codon ##

        ref_start = ref_dict[seq_id][gene_name][0]
        start_pos = -1
        end_pos = -1
        start_codon_missing = 1
        #if ref_start > 10:
        j_stop = ref_start + 10
        #else:
        #    j_stop = 10
        j = 0
        while j <= j_stop:
            sub_start_seq = ""
            if dire == '+':
                if (scaf_start-1) - j*3 >= 0:
                    start_pos = (scaf_start-1) - j*3
                    end_pos = (scaf_start-1) + 21
                    sub_start_seq = record.seq[start_pos:end_pos].upper()
                else:
                    break
            if dire == '-':
                if (scaf_end-1+1) + j*3 < seq_len:
                    start_pos = (scaf_end-1+1) - 21
                    end_pos = (scaf_end-1+1) + j*3
                    sub_start_seq = record.seq[start_pos:end_pos].upper()
                    sub_start_seq = sub_start_seq.reverse_complement()
                else:
                    break

            sub_start_seq += 'TAG'

            try:
                pro_seq_start = sub_start_seq.translate(table=table, cds=True, to_stop=False)
            except:
                pro_seq_start = sub_start_seq.translate(table=table, cds=False, to_stop=False)

            stop_count = pro_seq_start.count('*')
            #print(stop_count)
            if pro_seq_start.startswith('M') and stop_count <= 1:
                if dire == '+':
                    pos_revised_dict[seq_id][gene_name][0] = start_pos + 1
                    #print(start_pos + 1)
                elif dire == '-':
                    pos_revised_dict[seq_id][gene_name][1] = end_pos
                    #print(end_pos)
                #print(pro_seq_start)
                start_codon_missing = 0
                break

            j += 1


        ## if start codon missing, then try to find 'M' on backword direction.
        if start_codon_missing:
            j = 0
            while j < 10:
                if dire == '+':
                    start_pos =  scaf_start-1 + j*3
                    end_pos = scaf_start-1 + j * 3 + 3
                    sub_start_seq = record.seq[start_pos:end_pos].upper()
                if dire == '-':
                    start_pos = (scaf_end-1) -3 - j * 3 + 1
                    end_pos = (scaf_end-1) - j*3 + 1
                    sub_start_seq = record.seq[start_pos:end_pos].upper()
                    sub_start_seq = sub_start_seq.reverse_complement()

                sub_start_seq = sub_start_seq + 'TAG'
                try:
                    pro_seq_start = sub_start_seq.translate(table=table, cds=True, to_stop=False)
                except:
                    pro_seq_start = sub_start_seq.translate(table=table, cds=False, to_stop=False)

                if pro_seq_start.startswith('M'):
                    if dire == '+':
                        pos_revised_dict[seq_id][gene_name][0] = start_pos+1
                        #print(start_pos+1)
                    elif dire == '-':
                        #pos_revised_dict[seq_id][gene_name][0] = end_pos
                        pos_revised_dict[seq_id][gene_name][1] = end_pos
                        #print(end_pos)
                    #print(pro_seq_start)
                    start_codon_missing = 0
                    break

                j += 1

        #print(gene_name, start_codon_missing)
        ## if start codon still missing
        if start_codon_missing:
            start_still_missing_dict.setdefault(seq_id, []).append(gene_name)
            if dire == '+':
                #pos_revised_dict[seq_id][gene_name][0] = "<" + str(start_pos)
                pos_revised_dict[seq_id][gene_name][0] = "<" + str(scaf_start)
                pos_revised_dict[seq_id][gene_name][3] = "STARTNOTFOUND"
                # scaf_start, scaf_end
            elif dire == '-':
                #pos_revised_dict[seq_id][gene_name][1] = ">" + str(end_pos)
                pos_revised_dict[seq_id][gene_name][1] = ">" + str(scaf_end)
                pos_revised_dict[seq_id][gene_name][3] = "STARTNOTFOUND"

        ## find start codon ##
        ######################

###################################################################################################
STARTNOTFOUND = "can not determine the start codon"
STOPNOTFOUND = "can not determine the stop codon"

fh_out = open(out_f, 'w')
for seq_id in pos_revised_dict.keys():
    print(">%s" % seq_id, file=fh_out)
    pos_sort = {}
    for gene_name in pos_revised_dict[seq_id].keys():
        start, end, dire, TAA_note = pos_revised_dict[seq_id][gene_name]
        ## if start or end have '>' or '<' char, this may cause errors!
        if isinstance(start, str) or isinstance(end, str):
            if TAA_note == 'NA':
                print(gene_name, start, end, dire, sep="\t", file=fh_out)
            elif 'STARTNOTFOUND' in TAA_note:
                print(gene_name, start, end, dire, STARTNOTFOUND,
                     sep="\t", file=fh_out)
            elif 'STOPNOTFOUND' in TAA_note:
                print(gene_name, start, end, dire, STOPNOTFOUND,
                     sep="\t", file=fh_out)
            else:
                print(gene_name, start, end, dire, TAA_note, sep="\t", file=fh_out)
        else:
            pos_sort[start] = (gene_name, start, end, dire, TAA_note)
        #print(start)

    for key in sorted(pos_sort.keys()):
        #if sys.version_info > (3.0):
        print(pos_sort[key][0], pos_sort[key][1], pos_sort[key][2], pos_sort[key][3], sep="\t", end="\t", file=fh_out)
        if pos_sort[key][4] != 'NA':
            print(pos_sort[key][4], file=fh_out)
        else:
            print(file=fh_out)
        #else:
         #   print (pos_sort[key][0], "\t", pos_sort[key][1], "t\", pos_sort[key][2], "\t", pos_sort[key][3]),
          #  if pos_sort[key][4] != 'NA':
           #     print pos_sort[key][4]


if start_still_missing_dict:
    print("Following start codon still missing:")
    for seq_id in start_still_missing_dict.keys():
        print(">%s" % seq_id)
        for gene_name in start_still_missing_dict[seq_id]:
            print(gene_name)

if stop_still_missing_dict:
    print("Following stop codon still missing:")
    for seq_id in stop_still_missing_dict.keys():
        print(">%s" % seq_id)
        for gene_name in stop_still_missing_dict[seq_id]:
            print(gene_name)
