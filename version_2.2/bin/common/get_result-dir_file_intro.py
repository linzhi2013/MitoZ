#!/usr/bin/env python3
import sys
import os
from glob import glob
import re
from Bio import SeqIO


def get_file_intro(result_dir=None):

    result_dir = re.sub(r'/$', '', result_dir)
    files = glob(result_dir+'/*')

    kmer_mitofiles = []
    kmer_most_related_sp_files = []

    annt_most_related_sp_file = ''

    tbl2asn_gbf = ''
    tbl2asn_tbl = ''
    tbl2asn_fsa = ''
    tbl2asn_sqn = ''
    tbl2asn_val = []

    circos_image = []

    # gene sequences
    cds_file = ''
    rrna_file = ''
    whole_mt_file = ''
    trna_file = ''
    summary_file = ''

    contamination_seq_files = []

    for f in files:
        f = os.path.basename(f)
        if re.search(r'work\d+\.hmmout\.fa', f) or re.search(r'.*\.mitogenome\.fa', f):
            kmer_mitofiles.append(f)

        elif re.search(r'work.*most\_related\_species\.txt', f):
            kmer_most_related_sp_files.append(f)
            continue

        # must be after `kmer_most_related_sp_files`
        elif re.search(r'.*most\_related\_species\.txt', f):
            annt_most_related_sp_file = f

        elif re.search(r'\_mitoscaf\.fa\.tbl$', f):
            tbl2asn_tbl = f

        elif re.search(r'\_mitoscaf\.fa\.fsa', f):
            tbl2asn_fsa = f

        elif re.search(r'\_mitoscaf\.fa\.sqn', f):
            tbl2asn_sqn = f

        elif re.search(r'\_mitoscaf\.fa\.gbf', f):
            tbl2asn_gbf = f

        elif re.search(r'\_mitoscaf\.fa\.val', f):
            tbl2asn_val.append(f)

        elif re.search(r'errorsummary\.val', f):
            tbl2asn_val.append(f)

        elif re.search(r'circos\.(png|svg)', f):
            circos_image.append(f)

        elif re.search(r'\.trna', f):
            trna_file = f
        elif re.search(r'\.rrna', f):
            rrna_file = f
        elif re.search(r'\.cds', f):
            cds_file = f
        elif re.search(r'\.fasta', f) and not re.search(r'\_abundance', f):
            whole_mt_file = f
        elif re.search(r'summary.txt', f):
            summary_file = f

        elif re.search(r'\_abundance', f):
            contamination_seq_files.append(f)

    out_item = 0
    output = ''


    if (summary_file or tbl2asn_val or annt_most_related_sp_file):
        # add dash line
        output += '-' * 90
        output += '\n'
        output += '-' * 25
        output += ' summary information '
        output += '-' * 44
        output += '\n'

    if summary_file:
        out_item += 1
        output += '''
- {out_item}. The summary information for the whole analysis.
File:
{0}
'''.format(summary_file, out_item=out_item)

    if len(tbl2asn_val) > 0:
        out_item += 1
        output += '''
- {out_item}. You should check the following files to confirm if there are potential
annotation and assembly problems, e.g., a InternalStop for protein coding
genes can be caused by assemlby (e.g., there are Ns), annotation, mutation,
or the sequence simply comes from nuclear mitochondrial DNA segments (NUMTs).
You may inspect its sequencing depth from the visualization result files.
File:
{0}
'''.format('\n'.join(tbl2asn_val), out_item=out_item)

    if annt_most_related_sp_file:
        out_item += 1
        output += '''
- {out_item}. The most closely related species for the mitochondrial genome:
File:
{0}
        '''.format(annt_most_related_sp_file, out_item=out_item)



    if (tbl2asn_gbf or whole_mt_file or len(circos_image)>0 or tbl2asn_sqn or tbl2asn_tbl):
        # add dash line
        output += '\n'
        output += '-' * 90
        output += '\n'
        output += '-' * 25
        output += ' mitochondrial genome files '
        output += '-' * 37
        output += '\n'


    if tbl2asn_gbf:
        out_item += 1
        output += '''
- {out_item}. The mitochondrial genome in GenBank format.
File:
{0}
'''.format(tbl2asn_gbf, out_item=out_item)

    if whole_mt_file:
        out_item += 1
        output += '''
- {out_item}. The mitochondrial genome in Fasta format.
File:
{0}
'''.format(whole_mt_file, out_item=out_item)

    if len(circos_image) > 0:
        out_item += 1
        output += '''
- {out_item}. The visualization of the mitochondrial genome.
File:
{0}
'''.format('\n'.join(circos_image), out_item=out_item)

    if tbl2asn_sqn:
        out_item += 1
        output += '''
- {out_item}. The mitochondrial genome in sqn format.
This file can be used to upload to GenBank database after you confirm
the file {0} has no annotation problem.
File:
{1}
'''.format(tbl2asn_gbf, tbl2asn_sqn, out_item=out_item)

    if tbl2asn_tbl:
        out_item += 1
        output += '''
- {out_item}. The feature table for the mitochondrial genome.
File:
{0}
'''.format(tbl2asn_tbl, out_item=out_item)



    if (cds_file or rrna_file or trna_file):
        # add dash line
        output += '\n'
        output += '-' * 90
        output += '\n'
        output += '-' * 25
        output += ' each gene sequences '
        output += '-' * 44
        output += '\n'

    # cds_file = ''
    #rrna_file = ''
    #whole_mt_file = ''
    # trna_file =
    if cds_file:
        out_item += 1
        output += '''
- {out_item}. The CDS sequences of each protein coding gene.
File:
{0}
'''.format(cds_file, out_item=out_item)

    if rrna_file:
        out_item += 1
        output += '''
- {out_item}. The rRNA gene sequences for each rRNA gene.
File:
{0}
'''.format(rrna_file, out_item=out_item)

    if trna_file:
        out_item += 1
        output += '''
- {out_item}. The tRNA gene sequences for each tRNA gene.
File:
{0}
'''.format(trna_file, out_item=out_item)



    if (len(kmer_mitofiles) > 0 or len(kmer_most_related_sp_files) > 0):
        # add dash line
        output += '\n'
        output += '-' * 90
        output += '\n'
        output += '-' * 25
        output += ' mitochondrial sequences from each kmer assembly '
        output += '-' * 16
        output += '\n'

    if len(kmer_mitofiles) > 0:
        out_item += 1
        output += '''
- {out_item}. The mitochondrial sequences from each Kmer assembly.
File:
{0}
        '''.format('\n'.join(kmer_mitofiles), out_item=out_item)

    if len(kmer_most_related_sp_files) > 0:
        out_item += 1
        output += '''
- {out_item}. The most closely related species of each sequence of each Kmer assembly.
File:
{0}
        '''.format('\n'.join(kmer_most_related_sp_files), out_item=out_item)



    if len(contamination_seq_files) > 0:
        # add dash line
        output += '\n'
        output += '-' * 90
        output += '\n'
        output += '-' * 5
        output += ' potential contamination or nuclear mitochondrial DNA segments (NUMTs) (if any) '
        output += '-' * 5
        output += '\n'

    if len(contamination_seq_files) > 0:
        out_item += 1
        output += '''
- {out_item}. The fasta files of potential contamination or numts,
and their gene information (from HMMER search, which are just rough results).
File:
{0}
        '''.format('\n'.join(sorted(contamination_seq_files)), out_item=out_item)


    # add dash line

    output += '\n'
    output += '-' * 37
    output += ' THE END '
    output += '-' * 43
    output += '\n'

    print(output)


def main():
    usage = '''
python3 {0} <result_dir>
'''.format(sys.argv[0])
    if len(sys.argv) != 2:
        sys.exit(usage)

    get_file_intro(result_dir=sys.argv[1])




if __name__ == '__main__':
    main()









