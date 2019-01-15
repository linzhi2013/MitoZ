#!/usr/bin/python3
"""
find_final_MITObim_assembly.py

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
import os

usage = """
Description
	To find out the final MITObim.pl's assembly file.
	Then copy the file to current directory
Usage
	python3 {0} <abspath_to_search> <sample name> <ref name> <linkfile_name>
""".format(sys.argv[0])

if len(sys.argv) != 5:
	print(usage)
	sys.exit(0)

path, sample, ref, linkfile_name = sys.argv[1:5]

allfiles = os.listdir(path)
maxdir_dict = {}
for line in allfiles:
	if line.startswith("iteration"):
		dir_num = line.replace("iteration", "")
		maxdir_dict[dir_num] = line

dir_num_list = [int(i) for i in maxdir_dict.keys()]
maxdir = max(dir_num_list)
maxdir = str(maxdir)

prefix = sample + '-' + ref
# testpool-Trans-reslut_out_testpool.unpadded.fasta
target_path = os.path.join(path, maxdir_dict[maxdir],
	prefix+'_assembly',
    prefix+'_d_results', prefix+'_out_'+ "AllStrains.unpadded.fasta")
	#prefix+'_d_results', prefix+'_out_'+sample+'.unpadded.fasta')

if os.path.exists(target_path):
	#os.system("ln -s %s %s " % (target_path, linkfile_name))
	os.system("cp %s %s " % (target_path, linkfile_name))
else:
	print("can not find %s!" % target_path)
