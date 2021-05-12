#!/usr/bin/env python
import sys
import os
import re


def read_fastaLike2(file, seqid_pattern='^>', maxrecords=-1):
    '''
    Every time return one record using the 'yield' function,
    the return record is a list, containing the 'seqid' line,
    and 'sequence' lines.

    By default, the 'seqid' line has regular pattern '^>',while the 'sequence'
    lines don't. Change this behavior with 'seqid_pattern' option.

    Usage:
    >>>records = read_fastaLike2('myFastaLikefile')
    >>>for rec in records:
    >>>    print('seqid:', rec[0])
    >>>    print('first seq line:', rec[1])

    '''

    with open(file, 'r') as fh:
        firstline = True
        count = 0
        rec = []
        for i in fh:
            i = i.rstrip()
            if re.search(seqid_pattern, i):
                if firstline:
                    firstline = False
                else:
                    yield rec

                count += 1
                if (maxrecords > 0) and (count > maxrecords):
                    raise StopIteration

                rec = []
                rec.append(i)

            else:
                rec.append(i)

        yield rec

        raise StopIteration


def get_gene_coor(file):
    '''
    return:
    A list:

    seqid_geneCoor
    '''

    seqid_geneCoor = {}

    records = read_fastaLike2(file, seqid_pattern='^>Feature')
    for rec in records:
        tmp = False
        for index, i in enumerate(rec):
            line = i.split('\t')

            if index == 0:
                seqid = i.split()[1]
                seqid_geneCoor.setdefault(seqid, [])
                continue

            if len(line) == 3 and line[2] == 'gene':
                start, end = line[0:2]
                tmp = True
                continue

            # the next line
            if tmp:
                geneName = line[-1]
                seqid_geneCoor[seqid].append((start, end, geneName))
                tmp = False

    return seqid_geneCoor





















