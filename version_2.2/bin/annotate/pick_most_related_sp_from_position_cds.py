#!/usr/bin/env python3
import sys
import re
import os
import collections
from ete3 import NCBITaxa


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
            if seqid not in seqid_sp_count:
                seqid_mostRelatedSp[seqid] = 'Not found'
                print(seqid, 'not in input cds.taxa file!', file=sys.stderr)
                continue

            sp_count = seqid_sp_count[seqid]
            sp_count_sorted = sorted(
                                sp_count.items(),
                                key=lambda x:x[1],
                                reverse=True)

            seqid_mostRelatedSp[seqid] = sp_count_sorted[0][0]

    return seqid_mostRelatedSp

def get_rank_dict(taxa_name=None):
    ncbi = NCBITaxa()
    name_dict = ncbi.get_name_translator([taxa_name])
    if not name_dict:
        ## try only the first word (which may be a genus name?)
        print("can not find taxid for", taxa_name, file=sys.stderr)
        taxa_name = taxa_name.split()
        if len(taxa_name) > 1:
            taxa_name = taxa_name[0]
            print("try to search %s instead..." % taxa_name, file=sys.stderr)
            name_dict = ncbi.get_name_translator([taxa_name])

        if not name_dict:
            print("can not find taxid for %s, maybe it's a misspelling.\n" % taxa_name , file=sys.stderr)
            return None

    lineage_taxid_list = ncbi.get_lineage(name_dict[taxa_name][0])

    rank_dict = dict()
    for rank in ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']:
        rank_dict[rank] = 'NA'

    for j in lineage_taxid_list:
        rank = ncbi.get_rank([j])[j]
        taxa = ncbi.get_taxid_translator([j])[j]
        if rank == 'kingdom':
            rank_dict['kingdom'] = taxa

        elif rank == 'phylum':
            rank_dict['phylum'] = taxa

        elif rank == 'class':
            rank_dict['class'] = taxa

        elif rank == 'order':
            rank_dict['order'] = taxa

        elif rank == 'family':
            rank_dict['family'] = taxa

        elif rank == 'genus':
            rank_dict['genus'] = taxa

        elif rank == 'species':
            rank_dict['species'] = taxa

        else:
            pass

    return rank_dict


def get_lineage_for_cds(cdspos=None):
    ranks = ['kingdom', 'phylum', 'class', 'order', 'family','genus','species']

    rank_count = dict()
    for rank in ranks:
        rank_count[rank] = 0

    name_list = []

    # >gi_NC_011755_ATP8_Nezara_viridula_52_aa_Ref1:52aa
    # >gi_NC_KT345703_COX2_Capricornis_thar_227_aa-D4_Ref75:112aa   [mRNA]  locus=C3191348:60:146:-
    #pat = re.compile(r'>gi_NC_\d+_\w+?_(.+)_\d+_aa.*_Ref\d+:\d+aa.*locus=(.+?):')
    pat = re.compile(r'>gi_.+?_.+?_\w+?_(.+)_\d+_aa.*_Ref\d+:\d+aa.*locus=(.+?):')

    cdspos_taxa_log = cdspos+".taxa"
    with open(cdspos, 'r') as fh, open(cdspos_taxa_log, 'w') as fh_taxalog:
        for i in fh:
            if not i.startswith(">"):
                continue

            match = re.match(pat, i)
            if match:
                #species = match.group(1).replace("_", " ").replace(" sp.", "")
                species = match.group(1).replace("_", " ")
                if species.endswith('.'):
                    species = species.split()[0]

                print(match.group(2), match.group(1), sep="\t", end="\t", file=fh_taxalog)

                # assuming the rank_dict['species'] can not be 'NA',
                # and 'species' variable is a real species level name
                rank_dict = get_rank_dict(taxa_name=species)
                if not rank_dict:
                    print(file=fh_taxalog)
                    continue

                for rank in ranks:
                    print(rank_dict[rank], sep="|", end="|", file=fh_taxalog)
                    ## because some level of rank_dict may be 'NA'

                print(file=fh_taxalog)


            else:
                print("regular match fails:", i, flush=True, file=sys.stderr)
                continue

    return cdspos_taxa_log


def main():
    usage = '''
python3 {0} <work71.hmmtblout.besthit.sim.fa.solar.genewise.gff.cds.position.cds> <selected.fa> <outfile>
    '''.format(sys.argv[0])

    if len(sys.argv) != 4:
        sys.exit(usage)

    cds_pos_file, fa_file, outfile = sys.argv[1:4]

    cds_taxa_file = get_lineage_for_cds(cdspos=cds_pos_file)

    seqid_sp_count = get_seqid_taxa(cds_taxa_file=cds_taxa_file)

    seqid_mostRelatedSp = find_most_related_sp_for_each_seq(fa_file=fa_file,
        seqid_sp_count=seqid_sp_count)

    with open(outfile, 'w') as fhout:
        for seqid, sp in seqid_mostRelatedSp.items():
            print(seqid, sp, sep='\t', file=fhout)



if __name__ == '__main__':
    main()
