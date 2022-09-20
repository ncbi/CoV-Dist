#!/usr/bin/env python
# -*- coding: utf-8 -*-

import re
import sys
from collections import defaultdict
import csv

path_to_tsv=sys.argv[1]
out_path=sys.argv[2]

#path_to_tsv="20220918/lineagePaths.txt"
#out_path="ref_vcf"

lineageInfo=defaultdict()

def mutList2dt(lineage,mutList):
    dt=defaultdict()
    for line in mutList:
        m = re.match(r"(?P<ori>.)(?P<pos>\d+)(?P<der>.)", line)
        pos=m.group('pos')
        ori=m.group('ori')
        der=m.group('der')
        if pos in dt:
            if dt[pos][1] == ori:
                dt[pos][1] = der
            else:
                print("ERROR: pos:{}".format(pos))
            if dt[pos][0] == dt[pos][1]:
                del dt[pos]
        else:
            ele={pos:[ori,der]}
            dt.update(ele)
    return dt

def line2mutList(line):
    return [i for i in line.replace("> ",",").split(",") if i]


with open(path_to_tsv, mode='r') as csv_file:
    reader = csv.reader(csv_file, delimiter = '\t',skipinitialspace=True)
    for row in reader:
        lineage = row[0]
        if lineage == "19A" or lineage == "B" or lineage=="clade":
            pass
        else:
            mutList = line2mutList(row[2])
            dt = mutList2dt(lineage,mutList)
            lineageInfo.update({lineage:dt})

for i in lineageInfo:
    f=open(out_path+"/"+i+".vcf","w")
    f.write("#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT\n")
    for j in lineageInfo[i]:
        f.write("NC_045512.2\t"+str(j)+"\t.\t"+lineageInfo[i][j][0]+"\t"+lineageInfo[i][j][1]+"\t.\tPASS\t.\t.\n")
    f.close()
