#!/usr/bin/evn python3
import sys
import os
import re
import copy
import argparse
import collections
from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline
import subprocess

'''
This script is to search for missing PCGs from multiKmer assemlby results.

The script accepts the HMMER result sequences of multiKmer assemlbies, along
with their PCG annotation information output by HMMER search.

The script then tries to find the missing PCG with longest annotation length
from the HMMER result sequences, assuming the longer gene can be more likely
to be true positive genes.


If the new derived sequences also contain PCGs already present on the quick
mode source mito sequences, the script then blastn the new derived sequences
against the quick mode source mito sequences. The new derived sequences will
be further removed if it has low similarity (or even no match results) with
the quick mode source mito sequences, assuming such kind of new derived
sequences can be NUMTs or contaminations.

The new derived sequences which do not contain PCGs already present on the
quick mode source mito sequences will be kept.

-----------------------------------------------------------------------------
The clean fasta sequences output by this script still need to do PCG annotation
, only the sequence annotated with corresponding PCGs will be kept. If a
sequence also is also annotated with PCGs already present on the quick mode
sequence, and they have overlapping regions, Blastn this new drived sequence
to corresponding quick mode sequence,

For the sequences remained with PCGs, we should further check InternalStop,
alignment with sanger, similarity with sanger (and NCBI), and get CIGARS.

These precedures can further delete false positive results.


Copyright 2018. Guanliang MENG

'''


def get_parameter():
    desc = ''' To pick out the missing PCGs from multiKmer results '''

    parser = argparse.ArgumentParser(description=desc)

    parser.add_argument('--quick_mode_seq_file', metavar='<str>',
        required=False,help='quick mode fasta file')

    parser.add_argument('--quick_mode_fa_genes_file', metavar='<str>',
        required=False, help='gene name for each fasta sequence. format: seqid\t gene1 gene2')

    parser.add_argument('--multiKmer_sorted_scores_files', metavar='<str>',
                        nargs='+', help='files: work*.hmmtblout.besthit.sim.filtered.reformat.sorted')

    parser.add_argument('--multiKmer_seqs_files', metavar='<str>', nargs='+',
                        help='fasta files: work*.hmmtblout.besthit.sim.filtered.fa')

    parser.add_argument('--missing_genes', metavar='<str>', nargs='+',
                        help='missing gene names')

    parser.add_argument('--outprefix', metavar='<str>', required=True,
                        help='output fasta file name')

    parser.add_argument('--blastn', metavar='<str>', default='blastn',
                        help='path to the executable of blastn [%(default)s]')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()

    args = parser.parse_args()

    return args


def read_fastaLike(file=None):
    '''
    Every time return one record using the 'yield' function,
    which is a list, containing the 'seqid' line,
    and 'sequence' lines.

    Usage:
    >>>records = read_fastaLike('myFastaLikefile')
    >>>for rec in records:
    >>>    print('seqid:', rec[0])
    >>>    print('first seq line:', rec[1])

    '''

    with open(file, 'r') as fh:
        firstline = True
        rec = []
        for i in fh:
            i = i.rstrip()
            if i.startswith('>'):
                if firstline:
                    firstline = False
                else:
                    yield rec

                rec = []
                rec.append(i)

            else:
                rec.append(i)

        yield rec

        raise StopIteration


def get_hmmer_gene_annoations(multiKmer_sorted_scores_files=None):
    '''
    >scaffold512 Locus_1222_0 8.3 LINEAR length=1717 score=20.785
    COX2    2   649 45  643 +   4
    COX3    897 1691    18  784 +   4
    '''
    hmm_gene_annotations = collections.OrderedDict()
    for f in multiKmer_sorted_scores_files:
        # file name format:
        # work91.hmmtblout.besthit.sim.filtered.reformat.sorted
        # for fasta file:
        # work91.hmmtblout.besthit.sim.filtered.fa
        prefix = os.path.basename(f).split('.')[0]
        for rec in read_fastaLike(f):
            seqid = rec[0].split()[0].replace('>', '')
            if prefix not in hmm_gene_annotations:
                hmm_gene_annotations.setdefault(prefix, {})

            hmm_gene_annotations[prefix][seqid] = rec

    gene_and_seqids = collections.OrderedDict()
    for prefix in hmm_gene_annotations:
        for seqid in hmm_gene_annotations[prefix]:
            for index, item in enumerate(hmm_gene_annotations[prefix][seqid]):
                if index == 0:
                    continue
                gene, start, end = item.split()[0:3]
                start, end = [int(j) for j in (start, end)]
                gene_len = abs(end-start)
                if gene not in gene_and_seqids:
                    gene_and_seqids.setdefault(gene, [])

                gene_and_seqids[gene].append((gene, prefix, seqid, gene_len))

    return hmm_gene_annotations, gene_and_seqids


def pick_longest_seq_of_each_missing_gene(missing_genes=None, hmm_gene_annotations=None, gene_and_seqids=None, outprefix=None):
    res = collections.OrderedDict()

    gene_and_seqids_picked = []
    missing_genes = set(missing_genes)
    genes_found = set()
    for gene in gene_and_seqids:
        L = gene_and_seqids[gene]
        L_sorted = sorted(L, key=lambda x: x[3], reverse=True)
        gene, prefix, seqid, gene_len = L_sorted[0]
        if gene in missing_genes:
            genes_found.add(gene)
            gene_and_seqids_picked.append(L_sorted[0])

    genes_not_found = missing_genes - genes_found
    if len(gene_and_seqids_picked) > 0:
        print('\n{0} missing genes have been found from multiKmer scores files:'.format(
            len(gene_and_seqids_picked)), file=sys.stdout)
        for i in gene_and_seqids_picked:
            print('\t', i, sep='', file=sys.stdout)

        if len(genes_not_found) > 0:
            print('{0} missing genes still can not be found from multiKmer scores files:'.format(
                len(genes_not_found)), file=sys.stdout)
            print('\t', '\n'.join(genes_not_found), sep='', file=sys.stdout)
    else:
        raise ValueError(
            'No missing genes can be found from multiKmer scores files.')

    hmm_gene_annotations_picked = collections.OrderedDict()
    for i in gene_and_seqids_picked:
        prefix, seqid = i[1:3]
        if prefix not in hmm_gene_annotations_picked:
            hmm_gene_annotations_picked.setdefault(prefix, {})
        if seqid not in hmm_gene_annotations_picked[prefix]:
            hmm_gene_annotations_picked[prefix].setdefault(seqid, [])
            hmm_gene_annotations_picked[prefix][seqid].extend(
                hmm_gene_annotations[prefix][seqid])

    with open(outprefix+'.multiKmer.hmm.annotation', 'w') as fhout:
        for prefix in hmm_gene_annotations_picked:
            for seqid in hmm_gene_annotations_picked[prefix]:
                # add the prefix to the seqid line
                print(re.sub(r'^>', '>{0}_'.format(
                    prefix), hmm_gene_annotations_picked[prefix][seqid][0]), file=fhout)
                print(
                    '\n'.join(hmm_gene_annotations_picked[prefix][seqid][1:]), file=fhout)

    return hmm_gene_annotations_picked


def extract_multiKmer_seq_picked(multiKmer_seqs_files=None, hmm_gene_annotations_picked=None, seqformat='fasta', outfile=None):
    fhout = open(outfile, 'w')
    for prefix in hmm_gene_annotations_picked:
        for f in multiKmer_seqs_files:
            if os.path.basename(f).startswith(prefix):
                for rec in SeqIO.parse(f, seqformat):
                    seqid = str(rec.id)
                    # add the prefix to the seqid line
                    if seqid in hmm_gene_annotations_picked[prefix]:
                        seqid = '>{0}_{1}'.format(prefix, rec.description)
                        print(seqid, str(rec.seq), sep='\n', file=fhout)
                        # SeqIO.write(rec, fhout, seqformat)

                break
        else:
            print("not found multiKmer assembly file with prefix '{0}'".format(
                prefix), file=sys.stderr)

    fhout.close()

    return outfile


def determine_if_blastn(quick_mode_fa_genes_file=None, hmm_gene_annotations_picked=None):
    '''
    if at least one multiKmer sequence not only contain the missing genes,
    but also contain the genes already in quick mode results,
    then we need to do blastn
    '''

    if not quick_mode_fa_genes_file:
        # not providing quick mode results
        return False, None

    quick_mode_fa_genes = {}
    quick_mode_all_genes = []
    with open(quick_mode_fa_genes_file, 'r') as fh:
        for i in fh:
            i = i.strip()
            if not i:
                continue
            line = i.split()
            seqid = line[0]
            if seqid not in quick_mode_fa_genes:
                quick_mode_fa_genes.setdefault(seqid, [])
            for gene in line[1:]:
                quick_mode_all_genes.append(gene)
                quick_mode_fa_genes[seqid].append(gene)

    # hmm_gene_annotations[prefix][seqid] = rec
    # each 'rec' has content:
    # >work71_scaffold3477 Locus_8760_1 9.9 FORK length=725 score=12.146
    # ND6 7   126 338 423 +   2
    # ND2 114 725 20  570 +   4

    need_blastn = False
    gene_prefix_hmmseqid_quickmodeseqid_similarity_alnlen = []
    for prefix in hmm_gene_annotations_picked:
        for seqid in hmm_gene_annotations_picked[prefix]:
            rec = hmm_gene_annotations_picked[prefix][seqid]
            # the first line is seqid line, e.g. '>scaffold'
            for i in rec[1:]:
                gene = i.split()[0]
                for quickmodeseqid in quick_mode_fa_genes:
                    if gene in quick_mode_fa_genes[quickmodeseqid]:
                        need_blastn = True
                        tmp = [gene, prefix, seqid, quickmodeseqid, -1, -1]
                        gene_prefix_hmmseqid_quickmodeseqid_similarity_alnlen.append(
                            tmp)

    return need_blastn, gene_prefix_hmmseqid_quickmodeseqid_similarity_alnlen


def blastn_quickMode_and_multiKmer_seq_picked(quick_mode_seq=None, multiKmer_seq_picked=None, blastn=None, outfile='blastn.out'):

    makeblastdb = re.sub(r'blastn$', 'makeblastdb', blastn)
    cmd = makeblastdb + ' -in ' + quick_mode_seq + ' -dbtype nucl'
    subprocess.check_output(cmd, shell=True)

    blastn_cline = NcbiblastnCommandline(
        cmd=blastn,
        query=multiKmer_seq_picked,
        db=quick_mode_seq,
        out=outfile,
        outfmt=6,
        evalue=1e-5
    )

    blastn_cline()

    cmd = 'rm -rf {0}.{{{1}}}'.format(
        os.path.basename(quick_mode_seq), 'nin,nhr,nsq')
    subprocess.check_output(cmd, shell=True)

    return outfile


def get_lowsimilarity_query(blastn_tab_file=None, aln_len_cutoff=100, similarity_cutoff=98, gene_prefix_hmmseqid_quickmodeseqid_similarity_alnlen=None):
    '''
    Assumptions:
    multiKmer sequences with low similarities (or even no match) to quick mode
    sequences can be contaminations, numts.

    '''

    gene_prefix_hmmseqid_quickmodeseqid_similarity_alnlen_copy = copy.deepcopy(gene_prefix_hmmseqid_quickmodeseqid_similarity_alnlen)

    #print('Raw:', file=sys.stderr)
    #for j in gene_prefix_hmmseqid_quickmodeseqid_similarity_alnlen:
    #    print(j, file=sys.stderr)

    with open(blastn_tab_file, 'r') as fh:
        for i in fh:
            i = i.rstrip()
            query, sbjct, similarity, aln_len = i.split('\t')[0:4]
            similarity = float(similarity)
            aln_len = int(aln_len)

            prefix, seqid = query.split('_', maxsplit=1)
            for index, rec in enumerate(gene_prefix_hmmseqid_quickmodeseqid_similarity_alnlen_copy):
                # update the list with the top hit
                if (prefix in rec) and (seqid in rec) and (sbjct in rec) and gene_prefix_hmmseqid_quickmodeseqid_similarity_alnlen[index][-2] == -1:
                    gene_prefix_hmmseqid_quickmodeseqid_similarity_alnlen[index][-2] = similarity
                    gene_prefix_hmmseqid_quickmodeseqid_similarity_alnlen[index][-1] = aln_len
                else:
                    continue


    lowsimilarity_query = set()
    for rec in gene_prefix_hmmseqid_quickmodeseqid_similarity_alnlen:
        #print('\nupdated:')
        print(rec, file=sys.stderr)
        gene, prefix, hmmseqid, quickmodeseqid, similarity, aln_len = rec
        query = prefix+'_'+hmmseqid
        if (aln_len >= aln_len_cutoff) and (similarity < similarity_cutoff):
            # if has top-hit match meet the criteria
            lowsimilarity_query.add(query)
        elif (similarity == -1) and (aln_len == -1):
            # if no match, will also be treated as low similarity!
            lowsimilarity_query.add(query)

    if len(lowsimilarity_query) > 0:
        print('\nthe following sequences will be discarded:', file=sys.stderr)
        print('\n'.join(lowsimilarity_query), file=sys.stderr)
        return lowsimilarity_query

    else:
        print('\nNo low similarity multiKmer query', file=sys.stderr)

    return None


def filter_lowsimilarity_multiKmer(lowsimilarity_multiKmer_query=None, multiKmer_seq_picked_file=None, outfile=None):
    with open(multiKmer_seq_picked_file, 'r') as fh, open(outfile, 'w') as fhout:
        for rec in read_fastaLike(multiKmer_seq_picked_file):
            seqid_raw = re.sub('^>', '', rec[0].split()[0])
            if seqid_raw not in lowsimilarity_multiKmer_query:
                print('\n'.join(rec), file=fhout)

            #prefix, seqid = seqid_raw.split('_', maxsplit=1)

    return outfile


def main():
    args = get_parameter()

    hmm_gene_annotations, gene_and_seqids = get_hmmer_gene_annoations(
        multiKmer_sorted_scores_files=args.multiKmer_sorted_scores_files)

    hmm_gene_annotations_picked = pick_longest_seq_of_each_missing_gene(
        missing_genes=args.missing_genes,
        hmm_gene_annotations=hmm_gene_annotations,
        gene_and_seqids=gene_and_seqids,
        outprefix=args.outprefix
    )

    multiKmer_seq_picked_file = extract_multiKmer_seq_picked(
        multiKmer_seqs_files=args.multiKmer_seqs_files,
        hmm_gene_annotations_picked=hmm_gene_annotations_picked,
        seqformat='fasta',
        outfile=args.outprefix+'.multiKmer_seq_picked.fa'
    )

    need_blastn, gene_prefix_hmmseqid_quickmodeseqid_similarity_alnlen = determine_if_blastn(
        quick_mode_fa_genes_file=args.quick_mode_fa_genes_file, hmm_gene_annotations_picked=hmm_gene_annotations_picked)

    if need_blastn:
        blastn_tab_file = blastn_quickMode_and_multiKmer_seq_picked(
            quick_mode_seq=args.quick_mode_seq_file,
            multiKmer_seq_picked=multiKmer_seq_picked_file,
            blastn=args.blastn,
            outfile=args.outprefix+'.blastn.tab'
        )
        lowsimilarity_multiKmer_query = get_lowsimilarity_query(
            blastn_tab_file=blastn_tab_file,
            aln_len_cutoff=100,
            similarity_cutoff=98,
            gene_prefix_hmmseqid_quickmodeseqid_similarity_alnlen=gene_prefix_hmmseqid_quickmodeseqid_similarity_alnlen
        )

        if lowsimilarity_multiKmer_query:
            filter_lowsimilarity_multiKmer(
                lowsimilarity_multiKmer_query=lowsimilarity_multiKmer_query,
                multiKmer_seq_picked_file=multiKmer_seq_picked_file,
                outfile=args.outprefix+'.multiKmer_seq_picked.clean.fa'
            )

        else:
            cmd = 'cp {0} {1}'.format(
                multiKmer_seq_picked_file, args.outprefix+'.multiKmer_seq_picked.clean.fa')
            subprocess.check_output(cmd, shell=True)
    else:
        cmd = 'cp {0} {1}'.format(
            multiKmer_seq_picked_file, args.outprefix+'.multiKmer_seq_picked.clean.fa')
        subprocess.check_output(cmd, shell=True)


if __name__ == '__main__':
    main()
