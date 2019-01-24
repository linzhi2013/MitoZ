#!/usr/bin/python3
import os
import sys
import argparse
import subprocess
from Bio import SeqIO

usage = """
Description

    Reorder your mitochondrial genome sequence so as to be same with
    your reference. this is potenially needed becasue assember do not
    break mitogenome in "0" position for some reasons.

    Additionally, you can give a set of primer which is close to 0
    position of mitogenome, (e.g. F15: CACCCTATTAACCACTCACG for human)
    so I can use this information to adjust the sequence order.

    It supports three modes: (i). -f, using a reference guide; (ii).
    -p, using a human mitogenome primer set in a certain list(F15 or
    F361). However, it assumes that your primer is identical to
    sequence, if not, it will try find primer in ReverseComplementary
    sequence.

    If your mitogenome is circular, final reordered mitogenoem will be
    merely adjusted link orientation, if not, program will add 100 N
    between portions.


Usage

    python3 {0}  -f mito.fasta -r ref.fasta
    python3 {0}  -f mito.fasta -p F15

""".format(sys.argv[0])

parser = argparse.ArgumentParser(
    description=usage,
    formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('-f', metavar='FILE', type=str, required=True,
                    dest='mito', help="your mitogenome")

group = parser.add_mutually_exclusive_group()

group.add_argument('-r', metavar='FILE', type=str,
                    dest='ref', help="reference fasta")
group.add_argument('-p', metavar='STR', type=str, dest='primer',
                   choices=['F15', 'F361'],
                   help="primer set for [F15, F361], this is"
                   + " only for human (the primer sequences"
                   + " are included in the script)")

parser.add_argument('-m', metavar='INT', type=int,
                    dest='mismatch', default=2,
                    help="mismatch threshod for primer anchoring")

if len(sys.argv) == 1:
    parser.print_help()
    sys.exit()

args = parser.parse_args()

#primer set for anchoring "0" position
primers = {'F15':'CACCCTATTAACCACTCACG',
           'F361':'ACAAAGAACCCTAACACCAGC'}

# ------------------------------------------
def check_program_involed(cmd):
    '''
    check program involed whether is executable!
    '''
    result = (
        subprocess.call(
            "type %s" % cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        == 0
    )
    if result:
        return False
    else:
        return True

def read_fasta2dict(file):
    dict = {}
    topo = {}
    for rec in SeqIO.parse(file, 'fasta'):
        seq = str(rec.seq)
        id = rec.id.split()[0]
        dict[id] = seq
        if "circular" in rec.description:
            topo[id] = "circular"
        elif "linear" in rec.description:
            topo[id] = "linear"
        else:
            print("can not find topology=[circular|linear] in seqID")
            exit()

    return (dict, topo)

def distance(s1, s2):
    code= {
            '-' : 0,
            'A' : 2,
            'T' : 3,
            'C' : 5,
            'G' : 7,
            'R' : 14,
            'Y' : 15,
            'M' : 20,
            'K' : 21,
            'S' : 35,
            'W' : 24,
            'H' : 30,
            'B' : 105,
            'V' : 70,
            'D' : 42,
            'N' : 210
    }

    match = 0
    for i in range(len(s1)):
        a = s1[i]
        b = s2[i]
        if a == b and a != "-" :
            match += 1
        elif code[a] + code[b] > 14 and code[b] % code[a] == 0:
            match += 1
        else:
            False
    mismatch = len(s1) - match
    return mismatch

def find_primer_binding(primer, seq, mismatch):
    potential_site = []
    primerstr = primers[primer]
    primerloca = int(primer.replace("F", ""))
    plen = len(primerstr)
    qlen = len(seq)
    for i in range(qlen-plen):
        tmp = seq[i:i+plen]
        if distance(tmp, primerstr) <= mismatch:
            potential_site.append(i)
    if len(potential_site) >=2:
        print("more than 2 positions anchored!"
             + " please change a primer!")
        exit()
    elif len(potential_site) == 0:
        return False
    elif potential_site[0] >= primerloca:
        return potential_site[0]

def reverseComplement(s):
    complement = { 'A' : 'T', 'G' : 'C', 'C' : 'G', 'T' : 'A',
                  'a' : 'T', 'g' : 'C', 'c' : 'G', 't' : 'A'}
    t = ''
    for base in s:
        t = complement[base] + t
    return t
#------------------------------------------

# read mitogenome fasta
queryseq, querytopo = read_fasta2dict(args.mito)

out = open(os.path.basename(args.mito) + ".reorder", 'w')

# reference mode
if args.ref:
    # read reference fasta
    (refseq, none) = read_fasta2dict(args.ref)

    if len(refseq.keys()) > 1:
        print("can not accept more than one ref sequence")
        exit()
    else:
        Ref = list(refseq.values())[0].replace(" ", "")
        #print(">ref\n" + Ref)

    # blast
    if check_program_involed('makeblastdb') or check_program_involed('blastn'):
        print("can not find makeblastdb or blastn in your $PATH")
        exit()

    makeblastdb = "makeblastdb -in " +\
                args.ref +\
                " -dbtype nucl" +\
                " -parse_seqids" +\
                " -out ref"
    blastn_cmd = "blastn -query " +\
                args.mito +\
                " -out seq.blast" +\
                " -db ref" +\
                " -outfmt 6" +\
                " -evalue 1e-10" +\
                " -max_target_seqs 1" +\
                " -num_threads 2"

    subprocess.call(makeblastdb, shell=True)
    subprocess.call(blastn_cmd, shell=True)

    # prase reorder.blast

    blocks_dict = {}

    with open("seq.blast", 'r') as b:
        for i in b:
            tmp = i.strip().split()
            f = (tmp[8], tmp[9])
            m = (tmp[0], tmp[6], tmp[7])
            blocks_dict[f] = m
    # sort blocks by ref's position
    sorted_blocks = sorted(blocks_dict, key=lambda k: k[0])
    your_new_mito = ""
    for b in sorted_blocks:
        sr = b[0]
        er = b[1]
        tmp = blocks_dict[b]
        qid = tmp[0]
        sq = int(tmp[1]) - 1
        eq = int(tmp[2])
        your_new_mito += queryseq[qid][sq:eq]
        topo = "topology=" + querytopo[qid]

    out.write(">" + qid + " " + topo +"\n" + your_new_mito)
    os.system("rm ref.n*") # remove ref db

# primer using from defined list
elif args.primer:
    print("reorder your mitogenome accroding F15 primer set:"
         + "CACCCTATTAACCACTCACG")
    for i in queryseq.keys():
        if querytopo[i] == "circular":
            target = find_primer_binding(args.primer, queryseq[i], args.mismatch)
            if target:
                print("Primer " + args.primer + " found!")
                primerloca = int(args.primer.replace("F", ""))
                init = target - primerloca + 1
                your_new_mito = queryseq[i][init:] + queryseq[i][:init]
                out.write(">" + i + " topology=circular\n" + your_new_mito)
            else:
                print("Can not find primer "
                      + args.primer
                      + " in your mitogenome,"
                      + " now scanning it in reverseComplement fasta...")
                revcom = reverseComplement(queryseq[i])
                target = find_primer_binding(args.primer, revcom, args.mismatch)
                if target:
                    print("Primer " + args.primer + " found!")
                    primerloca = int(args.primer.replace("F", ""))
                    init = target - primerloca + 1
                    your_new_mito = revcom[init:] + revcom[:init]
                    out.write(">" + i + " topology=circular reverseComplement\n" + your_new_mito)
                else:
                    print("neithor in reverseComplement fasta")

        else:
            # this is for liner
            target = find_primer_binding(args.primer, queryseq[i], args.mismatch)
            if target:
                print("Primer " + args.primer + " found!")
                primerloca = int(args.primer.replace("F", ""))
                init = target - primerloca + 1
                your_new_mito = queryseq[i][init:] + "N" * 100 + queryseq[i][:init]
                out.write(">" + i + " topology=linear\n" + your_new_mito)
            else:
                print("can not find primer "
                      + args.primer
                      + " in your mitogenome,"
                      + " now scanning it in reverseComplement fasta...")
                revcom = reverseComplement(queryseq[i])
                target = find_primer_binding(args.primer, revcom, args.mismatch)

                if target:
                    print("Primer " + args.primer + " found!")
                    primerloca = int(args.primer.replace("F", ""))
                    init = target - primerloca + 1
                    your_new_mito = revcom[init:] + "N" * 100 + revcom[:init]
                    out.write(">" + i + " topology=linear reverseComplement\n" + your_new_mito)
                else:
                    print("can not anchor primer site")


