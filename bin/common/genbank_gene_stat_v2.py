#!/usr/bin/env python3
import sys
import os
import re
import collections

from Bio import SeqIO


mt_geneCount_std = {
    'l-rRNA': 1,
    's-rRNA': 1,

    'ATP6': 1,
    'ATP8': 1,
    'COX1': 1,
    'COX2': 1,
    'COX3': 1,
    'CYTB': 1,
    'ND1': 1,
    'ND2': 1,
    'ND3': 1,
    'ND4': 1,
    'ND4L': 1,
    'ND5': 1,
    'ND6': 1,

    'tRNA-Leu': 2,
    'tRNA-Val': 1,
    'tRNA-Phe': 1,
    'tRNA-Pro': 1,
    'tRNA-Thr': 1,
    'tRNA-Glu': 1,
    'tRNA-Leu': 2,
    'tRNA-Ser': 2,
    'tRNA-His': 1,
    'tRNA-Arg': 1,
    'tRNA-Gly': 1,
    'tRNA-Lys': 1,
    'tRNA-Asp': 1,
    'tRNA-Ser': 2,
    'tRNA-Tyr': 1,
    'tRNA-Cys': 1,
    'tRNA-Asn': 1,
    'tRNA-Ala': 1,
    'tRNA-Trp': 1,
    'tRNA-Met': 1,
    'tRNA-Gln': 1,
    'tRNA-Ile': 1,
}

def to_1_leftmost(pos=None):
    m = re.search(r'(\d+)', str(pos))
    pos_old = int(m.group(1))
    pos_new = pos_old + 1

    pos_old = str(pos_old)
    pos_new = str(pos_new)

    return re.sub(pos_old, pos_new, str(pos))


def gene_stat(gbfile=None):
    # start stop type gene product onScaffoldID
    # then report duplicate genes and missing genes
    gene_freq = {}

    gene_infor = []
    for rec in SeqIO.parse(gbfile, 'gb'):
        seqid = str(rec.id)
        for fea in rec.features:
            if fea.type not in ['CDS', 'rRNA', 'tRNA']:
                continue

            start = fea.location.start
            start = to_1_leftmost(start)

            end = fea.location.end
            end = to_1_leftmost(end)

            strand = fea.location.strand  # is a number
            if strand == 1:
                strand = '+'
            if strand == -1:
                strand = '-'
            gene = ""
            product = ""

            if 'gene' in fea.qualifiers:
                gene = fea.qualifiers['gene'][0]

            if 'product' in fea.qualifiers:
                product = fea.qualifiers['product'][0]
                if not gene:
                    gene = product
            else:
                print(ass_num, "Warning: NO gene or product tag! this gene is not output!\n")
                continue
            gene_for_count = gene
            if fea.type == 'tRNA':
                gene_for_count = product
            if gene_for_count not in gene_freq:
                gene_freq[gene_for_count] = 1
            else:
                gene_freq[gene_for_count] += 1
            line = [str(j) for j in (seqid, start, end, strand, fea.type, gene, product)]
            gene_infor.append(line)

    return gene_infor, gene_freq


def get_dup_and_missing_genes(gene_freq=None, MitoZ_module_used=None):
    PCGs_found = 0
    tRNA_found = 0
    rRNA_found = 0
    duplicted_genes = []
    for gene in gene_freq:
        count = gene_freq[gene]
        # duplication check
        if count > 1:
            duplicted_genes.append((gene, count))

        mt_geneCount_std[gene] -= count

        if not re.search(r'RNA', gene):
            # is PCG
            PCGs_found += count
        elif re.search(r'rRNA', gene):
            rRNA_found += count
        elif re.search(r'tRNA', gene):
            tRNA_found += count

    print('\n\n{:25s}{:>3}'.format('Protein coding genes totally found:', PCGs_found))
    print('{:33s}{:>5}'.format('tRNA genes totally found:', tRNA_found))
    print('{:33s}{:>5}'.format('rRNA genes totally found:', rRNA_found))
    print('-'*38)
    print('{:33s}{:>5}'.format('Genes totally found:', sum([rRNA_found,tRNA_found, PCGs_found])))

    missing_genes = []
    for gene in mt_geneCount_std:
        count = mt_geneCount_std[gene]
        if count > 0:
            missing_genes.append((gene, count))

    #if len(duplicted_genes) > 0 or len(missing_genes) > 0:
    #if len(missing_genes) > 0:
    #    print('\n\nWhen assuming there are 13 protein coding genes,\n22 tRNA genes (including 2 tRNA-Leu genes and two tRNA-Ser genes\nin many species), and 2 rRNA genes for one mitochondrial genome:')
    # if len(duplicted_genes) > 0:
    #    print('\nPotential duplicate genes:')
    #    print('{:20s}{:>15}'.format('#Gene', 'total_occur_number'))
    #    print('-'*38)
    #    for gene, count in sorted(duplicted_genes, key=lambda x:x[0]):
    #        print('{:33s}{:>5}'.format(gene, count))

    if len(missing_genes) > 0:
        print('\n\nPotential missing genes:')
        print('{:18s}{:>15}'.format('#Gene', 'total_missing_number'))
        print('-'*38)
        for gene, count in sorted(missing_genes, key=lambda x:x[0]):
            print('{:33s}{:>5}'.format(gene, count))

        if MitoZ_module_used in ['all', 'all2']:
            print("\nThe missing genes might be foud from the\n'*.high_abundance*' and '*.low_abundance*' files!")

        if MitoZ_module_used in ['annotate']:
            print("\nIf you got the mitochondrial genome using the `findmitoscaf`\nor `assemble` module of MitoZ, the missing genes might be\nfound from the '*.high_abundance*' and '*.low_abundance*' files\nin the 'tmp' directory!")




def get_seq_topology_and_related_sp(gbfile=None, closely_related_sp_file=None):
    seqid_len_topology_relatedSP = collections.defaultdict(list)
    for rec in SeqIO.parse(gbfile, 'gb'):
        seqlen = str(len(rec))
        seqid = str(rec.id)
        topology = rec.annotations['topology']
        circular_yes_no = 'no'
        if topology == 'circular':
            circular_yes_no = 'yes'

        seqid_len_topology_relatedSP[seqid].extend([seqlen, circular_yes_no])

    with open(closely_related_sp_file, 'r') as fh:
        for i in fh:
            i = i.strip()
            seqid, sp = i.split('\t')
            seqid_len_topology_relatedSP[seqid].append(sp)

    return seqid_len_topology_relatedSP


def main():
    usage = '''
python3 {0} <in.gb> <*.most_related_species.txt> <MitoZ_module_used>
    '''.format(sys.argv[0])

    if len(sys.argv) != 4:
        sys.exit(usage)

    gbfile, closely_related_sp_file, MitoZ_module_used = sys.argv[1:4]

    seqid_len_topology_relatedSP = get_seq_topology_and_related_sp(gbfile=gbfile, closely_related_sp_file=closely_related_sp_file)

    col_fomrat = '{:15s}{:15s}{:15s}{:20s}'
    tpo_title_line = col_fomrat.format('#Seq_id', 'Length(bp)', 'Circularity', 'Closely_related_species')
    print(tpo_title_line)

    for seqid in seqid_len_topology_relatedSP:
        seqlen, topology, sp = seqid_len_topology_relatedSP[seqid]
        out = col_fomrat.format(seqid, seqlen, topology, sp)
        print(out)

    print("\n")

    gene_infor, gene_freq = gene_stat(gbfile)
    col_fomrat = '{:15s}{:7s}{:7s}{:11s}{:11s}{:7}{:11s}{:34}{:15s}'
    title_line = col_fomrat.format('#Seq_id', 'Start','End', 'Length(bp)', 'Direction', 'Type','Gene_name','Gene_prodcut','Total_freq_occurred')
    print(title_line)
    print('-'*122)
    for i in gene_infor:
        start, end = i[1:3]
        start = int(re.sub(r'>|<', '', str(start)))
        end = int(re.sub(r'>|<', '', str(end)))
        gene_len = str(end-start+1)
        gene_type = i[4]
        gene_for_count = i[5]
        if gene_type == 'tRNA':
            gene_for_count = i[6]

        freq = gene_freq[gene_for_count]

        content_line = col_fomrat.format(i[0], i[1], i[2], gene_len, i[3], i[4], i[5], i[6], str(freq))
        print(content_line)

    get_dup_and_missing_genes(gene_freq=gene_freq, MitoZ_module_used=MitoZ_module_used)

if __name__ == '__main__':
    main()










