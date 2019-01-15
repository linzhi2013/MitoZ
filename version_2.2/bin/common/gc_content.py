#!/usr/bin/python
"""
gc_content.py

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
def gc_content(seq, window_size):
    seq = seq.upper()
    length = len(seq)
    L = []
    G_count = 0
    C_count = 0
    i = 0
    while i <= (length-window_size):
        subseq = seq[i:i+window_size]
       # print(i, subseq)
        G_count = subseq.count('G')
        C_count = subseq.count('C')
        GC_rate = (G_count+C_count)*100/window_size
        L.append((i+1,GC_rate))
        i += 1

    GC_total = 0
    GC_total = seq.count('G')
    GC_total += seq.count('C')
    GC_avg = GC_total*100/length

    return L, GC_avg

if __name__ == '__main__':
    print(gc_content(sys.argv[1], int(sys.argv[2])))
