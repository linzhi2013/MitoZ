### Name
usage: Mitogenome_reorder.py [-h] -f FILE [-r FILE | -p STR] [-m INT]

### Description

    Reorder your mitochondrial genome sequence so as to be same with
    your reference. this is potenially needed becasue assember do not
    break mitogenome in "0" position for some reasons.

    Additionally, you can give a set of primer which is close to 0
    position of mitogenome, (e.g. F15: CACCCTATTAACCACTCACG for human)
    so I can use this information to adjust the sequence order.

    It supports three modes: (i). -f, using a reference guide; (ii).
    -p, using a human mitogenome primer set in a certain list(F15 or
    F361). However, it assumes that your primer is identical to
    sequence, complementary or reverse primer not acceptable!

    If your mitogenome is circular, final reordered mitogenoem will be
    merely adjusted link orientation, if not, program will add 100 N
    between portions.

### Usage

    python3 Mitogenome_reorder.py  -f mito.fasta -r ref.fasta
    python3 Mitogenome_reorder.py  -f mito.fasta -p L15

### Optional arguments:

	  -h, --help  show this help message and exit
	  -f FILE     your mitogenome fasta
	  -r FILE     reference fasta
	  -p STR      primer set for [F15, F361], this is only for human
	  -m INT      mismatch threshod for primer anchoring

	  
### Tips

- Make sure you sequence ID contains "topology=circular" or "topology=linear", this is important!
- If you want to use your own sequnce as primer bait, you can modify this python script, add you priemr sequence as a pair of KEY-VALUSE to dict of "primers", and make sure the numer is accurate, e.g. F15 means that primer starts at 15th base in sequnce.
- Primer sequence can contain a degenerate base, like R, Y, M ...