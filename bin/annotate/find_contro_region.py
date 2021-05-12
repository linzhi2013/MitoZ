#!/usr/bin/env python3
import sys
import os
import re
import collections
from Bio import SeqIO
from itertools import groupby
import argparse

sys.path.insert(0, os.path.dirname(__file__))
from gene_feature import get_gene_coor


def get_para():
    desc = '''to find the control region'''

    parser = argparse.ArgumentParser(description=desc)

    parser.add_argument('-fa_file', metavar='file', help='input fasta file')

    parser.add_argument('-PCG_cutoff_file', metavar='file', help='PCG length cutoff file')

    parser.add_argument('-PCG_len_ratio', metavar='float', default=0.8,
        type=float,
        help='the PCG length must be larger than this ratio. [%(default)s]')

    parser.add_argument('-s_rRNA_CM_file', metavar='file',
        help='s-rRNA CM file')

    parser.add_argument('-l_rRNA_CM_file', metavar='file',
        help='l-rRNA CM file')

    parser.add_argument('-rRNA_len_ratio', metavar='float', default=0.8,
        type=float,
        help='the rRNA length must be larger than this ratio. [%(default)s]')

    parser.add_argument('-tRNA_num_min', metavar='int', default=22, type=int,
        help='the minimum number of tRNA on one seq [%(default)s]')

    parser.add_argument('-fea_files', metavar='file', nargs='+',
        help='feature table files')

    parser.add_argument('-CR_len_min', metavar='int', type=int, default=600,
        help='the minimum non-coding length to do control region annotation [%(default)s]')

    parser.add_argument('-outfile', metavar='file', help='output file')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    return parser.parse_args()


def read_rRNA_CM_length(CM_files=None, ratio=1):
    rRNA_CM_len = {}
    for CM_file in CM_files:
        if '12S' in CM_file.upper():
            gene = 's-rRNA'
        elif '16S' in CM_file.upper():
            gene = 'l-rRNA'

        with open(CM_file, 'r') as fh:
            for i in fh:
                i = i.strip()
                if not i:
                    continue
                m = re.match(r'^CLEN\s+(\d+)$', i)
                if m:
                    CM_len = m.group(1)
                    rRNA_CM_len[gene] = int(CM_len) * ratio
                    break

    return rRNA_CM_len


def read_PCG_length_cutoff(CDS_length_file=None, ratio=1):
    PCG_len = {}
    with open(CDS_length_file, 'r') as fh:
        for i in fh:
            i = i.strip()
            if not i:
                continue
            gene, length = i.split()[0:2]
            PCG_len[gene] = int(length)

    return PCG_len


def file_not_empty(file=None):
    """
    check if file is empty

    """
    if os.stat(file).st_size > 0:
        return True
    else:
        return False



def read_featureTable(fea_files=None):
    all_seqid_geneCoor = collections.defaultdict(list)

    for fea_f in fea_files:
        if not file_not_empty(fea_f):
            continue

        seqid_geneCoor = get_gene_coor(fea_f)
        # seqid_geneCoor.append((start, end, geneName))
        for seqid in seqid_geneCoor:
            for start, end, geneName in seqid_geneCoor[seqid]:
                all_seqid_geneCoor[seqid].append((start, end, geneName))

    return all_seqid_geneCoor


def confirm_geneCoor(all_seqid_geneCoor=None, PCG_len=None, rRNA_CM_len=None, tRNA_num=22):
    no_need_seqids = []
    for seqid in all_seqid_geneCoor:
        PCGs_met = set()
        tRNA_met = set()
        rRNA_met = set()


        for start, end, geneName in all_seqid_geneCoor[seqid]:
            start = re.sub(r'>|<', '', start)
            end = re.sub(r'>|<', '', end)

            start, end = [int(j) for j in (start, end)]
            if start > end:
                start, end = end, start
            gene_len = abs(end-start) + 1

            if geneName in PCG_len:
                # is PCG
                if PCG_len[geneName] <= gene_len:
                    PCGs_met.add(geneName)
                else:
                    print(seqid, geneName, 'too short: {0} bp < {1} bp'.format(gene_len, PCG_len[geneName]), file=sys.stderr)

            if geneName.startswith('trn'):
                # is tRNA gene
                tRNA_met.add(geneName)

            if geneName in rRNA_CM_len:
                # is rRNA
                if rRNA_CM_len[geneName] <= gene_len:
                    rRNA_met.add(geneName)
                else:
                    print(seqid, geneName, 'too short: {0} bp < {1} bp'.format(gene_len, rRNA_CM_len[geneName]), file=sys.stderr)

        if len(PCGs_met) < 13:
            print('Not enough 13 PCGs', file=sys.stderr)
            no_need_seqids.append(seqid)

        if len(tRNA_met) < tRNA_num:
            print('Not enough {0} tRNAs.'.format(tRNA_num), file=sys.stderr)
            no_need_seqids.append(seqid)

        if len(rRNA_met) < 2:
            print('Not enough 2 rRNAs.', file=sys.stderr)
            no_need_seqids.append(seqid)

    for seqid in no_need_seqids:
        all_seqid_geneCoor.pop(seqid, None)
        #del all_seqid_geneCoor[seqid]

    if len(all_seqid_geneCoor) > 0:
        return all_seqid_geneCoor
    else:
        return False


# only when the all_seqid_geneCoor has remains
# do the following
def get_seqid_lst(all_seqid_geneCoor=None, fa_file=None):
    seqid_lst = {}
    seqid_topology_and_len = {}

    for rec in SeqIO.parse(fa_file, 'fasta'):
        topology = 'linear'
        seqid = str(rec.id)
        if seqid not in all_seqid_geneCoor:
            continue
        seqid_lst[seqid] = [0] * len(rec)
        if 'topology=circular' in rec.description:
            topology = 'circular'
        seqid_topology_and_len[seqid] = (topology, len(rec))

    # all_seqid_geneCoor[seqid].append((start, end, geneName))
    for seqid in all_seqid_geneCoor:
        for start, end, geneName in all_seqid_geneCoor[seqid]:
            start = re.sub(r'>|<', '', start)
            end = re.sub(r'>|<', '', end)
            start, end = [int(j) for j in (start, end)]
            if start > end:
                start, end = end, start
            for x in range(start-1, end):
                seqid_lst[seqid][x] = 1
    return seqid_lst, seqid_topology_and_len


def group_regions(seqid_lst=None):
    seqid_start_len = collections.defaultdict(list)
    for seqid in seqid_lst:
        d = seqid_lst[seqid]
        start = 0
        # d = [1,1,1,0,0,1,1]
        # [[1, 1, 1], [0, 0], [1, 1]]
        for x in [list(g) for k, g in groupby(d)]:
            range_len = len(x)
            if x[0] == 0:
                seqid_start_len[seqid].append((start, range_len))
            start = start + range_len

    return seqid_start_len


def filter_range(seqid_start_len=None, range_len_cutoff=600, seqid_topology_and_len=None):
    seqid_start_len_filtered = collections.defaultdict(list)
    broken_CR = collections.defaultdict(list)

    for seqid in seqid_start_len:
        topology, seqlen = seqid_topology_and_len[seqid]
        x = seqid_start_len[seqid]
        # if it was circular
        # first check the start region and end region
        # which may be broken Control region
        if topology == 'circular':
            print('I am here')
            start1, range_len1 = x[0]
            start2, range_len2 = x[-1]
            if (range_len1 + range_len2) >= range_len_cutoff:
                print('join start + end ')
                broken_CR[seqid].append((start1,range_len1))
                broken_CR[seqid].append((start2,range_len2))

            for start, range_len in x[1:-1]:
                if range_len >= range_len_cutoff:
                    seqid_start_len_filtered[seqid].append((start, range_len))
        else:
            for start, range_len in x:
                if range_len >= range_len_cutoff:
                    seqid_start_len_filtered[seqid].append((start, range_len))

    return seqid_start_len_filtered, broken_CR


def main():

    args = get_para()

    all_seqid_geneCoor = read_featureTable(fea_files=args.fea_files)

    rRNA_CM_len = read_rRNA_CM_length(
        CM_files=[args.s_rRNA_CM_file, args.l_rRNA_CM_file],
        ratio=args.rRNA_len_ratio)

    PCG_len = read_PCG_length_cutoff(CDS_length_file=args.PCG_cutoff_file,
        ratio=args.PCG_len_ratio)

    # print(all_seqid_geneCoor)

    # filter some seqids
    all_seqid_geneCoor = confirm_geneCoor(
        all_seqid_geneCoor=all_seqid_geneCoor,
        PCG_len=PCG_len,
        rRNA_CM_len=rRNA_CM_len,
        tRNA_num=args.tRNA_num_min)

    # print(all_seqid_geneCoor)
    fhout = open(args.outfile, 'w')

    if not all_seqid_geneCoor:
        return None

    seqid_lst, seqid_topology_and_len = get_seqid_lst(
        all_seqid_geneCoor=all_seqid_geneCoor,
        fa_file=args.fa_file)

    # print(seqid_lst)
    if len(seqid_lst) == 0:
        return None

    # print(seqid_lst)
    seqid_start_len = group_regions(seqid_lst=seqid_lst)
    # print(seqid_start_len)

    seqid_start_len_filtered, broken_CR = filter_range(
        seqid_start_len=seqid_start_len,
        range_len_cutoff=args.CR_len_min,
        seqid_topology_and_len=seqid_topology_and_len)

    if len(seqid_start_len_filtered) > 0 and len(broken_CR) > 0:
        for seqid in seqid_start_len_filtered:
            # for each seqid, choose the longest non-coding
            # region as control region
            seqid_start_rangeLen = []
            for start, range_len in seqid_start_len_filtered[seqid]:
                seqid_start_rangeLen.append((seqid, start, range_len))

            seqid_start_rangeLen_sorted = sorted(seqid_start_rangeLen, key=lambda x: x[2], reverse=True)

            seqid, start, range_len = seqid_start_rangeLen_sorted[0]
            # if the sequence is circular
            # compare the two ends
            if len(broken_CR) > 0:
                start1, range_len1 = broken_CR[seqid][0]
                start2, range_len2 = broken_CR[seqid][1]
            else:
                range_len1 = range_len2 = 0

            if (range_len1 + range_len2) < range_len:
                begin = start + 1
                end = start + range_len

                print('>Feature {0}'.format(seqid), file=fhout)
                print(begin, end, 'misc_feature', sep='\t', file=fhout)
                print('\t\t\tnote', 'putative control region', sep='\t', file=fhout)
            else:
                print('>Feature {0}'.format(seqid), file=fhout)
                print(start1+1, start1+range_len1, 'misc_feature', sep='\t', file=fhout)
                print(start2+1, start2+range_len2, sep='\t', file=fhout)
                print('\t\t\tnote', 'putative control region', sep='\t', file=fhout)

    elif len(seqid_start_len_filtered) > 0 and len(broken_CR) == 0:
        for seqid in seqid_start_len_filtered:
            # for each seqid, choose the longest non-coding
            # region as control region
            seqid_start_rangeLen = []
            for start, range_len in seqid_start_len_filtered[seqid]:
                seqid_start_rangeLen.append((seqid, start, range_len))

            seqid_start_rangeLen_sorted = sorted(seqid_start_rangeLen, key=lambda x: x[2], reverse=True)

            seqid, start, range_len = seqid_start_rangeLen_sorted[0]

            begin = start + 1
            end = start + range_len

            print('>Feature {0}'.format(seqid), file=fhout)
            print(begin, end, 'misc_feature', sep='\t', file=fhout)
            print('\t\t\tnote', 'putative control region', sep='\t', file=fhout)

    elif len(seqid_start_len_filtered) == 0 and len(broken_CR) > 0:
        for seqid in broken_CR:
            start1, range_len1 = broken_CR[seqid][0]
            start2, range_len2 = broken_CR[seqid][1]
            print('>Feature {0}'.format(seqid), file=fhout)
            print(start1+1, start1+range_len1, 'misc_feature', sep='\t', file=fhout)
            print(start2+1, start2+range_len2, sep='\t', file=fhout)
            print('\t\t\tnote', 'putative control region', sep='\t', file=fhout)

    else:
        return None


    fhout.close()



if __name__ == '__main__':
    main()
