#!/usr/bin/env python3
import sys
import re
import os
import collections

'''
Query   Requiring_taxa  Metazoa|Chordata|NA|NA|NA|NA|NA|
C1      Pan_paniscus    Metazoa|Chordata|Mammalia|Primates|Hominidae|Pan|Pan paniscus|
C1      Homo_sapiens    Metazoa|Chordata|Mammalia|Primates|Hominidae|Homo|Homo sapiens|

'''

def get_seqid_taxa(cds_taxa_file=None):
    seqid_sp_count = collections.defaultdict(dict)
    with open(cds_taxa_file, 'r') as fh:
        for i in fh:
            i = i.strip()
            seqid, species = i.split('\t')[0:2]
            species = species.replace("_", " ")
            if species.endswith('.'):
                species = species.split()[0]

            if species not in seqid_sp_count[seqid]:
                seqid_sp_count[seqid][species] = 1
            else:
                seqid_sp_count[seqid][species] += 1

    return seqid_sp_count


def find_most_related_sp_for_each_seq(fa_file=None, seqid_sp_count=None):
    seqid_mostRelatedSp = collections.OrderedDict()
    with open(fa_file, 'r') as fh:
        for i in fh:
            i = i.strip()
            if not i.startswith('>'):
                continue
            seqid = i.split()[0].replace('>', '')
            if '_FivePCGs' in seqid:
                seqid_tmp = seqid.split('_FivePCGs')[0]
            else:
                seqid_tmp = seqid
            if seqid_tmp not in seqid_sp_count:
                seqid_mostRelatedSp[seqid] = 'Not found'
                print(seqid, 'not in input cds.taxa file!', file=sys.stderr)
                continue

            sp_count = seqid_sp_count[seqid_tmp]
            sp_count_sorted = sorted(
                                sp_count.items(),
                                key=lambda x:x[1],
                                reverse=True)

            seqid_mostRelatedSp[seqid] = sp_count_sorted[0][0]

    return seqid_mostRelatedSp


def main():
    usage = '''
python3 {0} <work71.hmmtblout.besthit.sim.fa.solar.genewise.gff.cds.position.cds.taxa> <selected.fa> <outfile>
    '''.format(sys.argv[0])

    if len(sys.argv) != 4:
        sys.exit(usage)

    cds_taxa_file, fa_file, outfile = sys.argv[1:4]

    seqid_sp_count = get_seqid_taxa(cds_taxa_file=cds_taxa_file)

    seqid_mostRelatedSp = find_most_related_sp_for_each_seq(fa_file=fa_file,
        seqid_sp_count=seqid_sp_count)

    with open(outfile, 'w') as fhout:
        for seqid, sp in seqid_mostRelatedSp.items():
            print(seqid, sp, sep='\t', file=fhout)



if __name__ == '__main__':
    main()
