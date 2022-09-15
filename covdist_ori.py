#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys, os, subprocess, shutil, csv
import argparse,  platform
import Bio.SeqIO, pysam, numpy as np
from time import time
from multiprocessing import Pool, cpu_count
from functools import partial
from itertools import combinations
from collections import OrderedDict

std_vec = {
        "A":[1.0,0.0,0.0,0.0],
        "C":[0.0,1.0,0.0,0.0],
        "G":[0.0,0.0,1.0,0.0],
        "T":[0.0,0.0,0.0,1.0],
        "a":[1.0,0.0,0.0,0.0],
        "c":[0.0,1.0,0.0,0.0],
        "g":[0.0,0.0,1.0,0.0],
        "t":[0.0,0.0,0.0,1.0]
        }

def check_path_existence(path):
    abspath = os.path.abspath(path)
    try:
        f = open(abspath)
    except IOError:
        sys.exit("File: {} not accessible! Please check!".format(abspath))
    finally:
        f.close()
    return(abspath)

def thetayc(v1, v2):
    # check if this position has coverage
    if np.sum(v1) == 0 or np.sum(v2) ==0:
        return 0
    else:
        # turn the vector into frequence 
        vec1 = np.true_divide(v1,np.sum(v1))
        vec2 = np.true_divide(v2,np.sum(v2))
        dotsum = np.dot(vec1,vec2)
        dist = np.sum(np.square(vec1),axis=0) + np.sum(np.square(vec2),axis=0) 
        if dist == dotsum:
            return 0
        else:
            return 1-dotsum/(dist-dotsum)

def round_vec(vec):
    return [round(i,3) for i in vec]

def CalculateDist(count_table1, count_table2):
    # with open("test_ori", "w") as g:
    #     for key in RefTableFreq:
    #         if key in count_table1:
    #             vec1 = round_vec(count_table1[key].count_vec)
    #         else:
    #             vec1 = round_vec(RefTableFreq[key].count_vec)
    #         if key in count_table2:
    #             vec2 = round_vec(count_table2[key].count_vec)
    #         else:
    #             vec2 = round_vec(RefTableFreq[key].count_vec)
    #         vec1 = np.true_divide(vec1,np.sum(vec1))
    #         vec2 = np.true_divide(vec2,np.sum(vec2))
    #         val = thetayc(vec1, vec2)
    #         v1 = "\t".join(list(map(str, vec1)))
    #         v2 = "\t".join(list(map(str, vec2)))
    #         g.write(str(key)+"\t"+v1+"\t"+v2+"\t"+str(val)+"\n")

    sumval = 0
    key_inter = set(count_table1).intersection(set(count_table2))
    for key in key_inter:
        sumval += thetayc(count_table1[key].count_vec, count_table2[key].count_vec)
    for key in (set(count_table1) - key_inter):
        sumval += thetayc(count_table1[key].count_vec, std_vec[count_table1[key].ref_allele])
    for key in (set(count_table2) - key_inter):
        sumval += thetayc(count_table2[key].count_vec, std_vec[count_table2[key].ref_allele])
    return sumval

class PosCount:
    def __init__(self,ref_id,ref_pos,ref_allele,count_vec):
        self.ref_id = ref_id
        self.ref_pos = ref_pos
        self.ref_allele = ref_allele
        self.count_vec = count_vec

    def __repr__(self):
        return "{0}(id={1!r}, pos={2!r}, count_vec={3!r})".format(
                self.__class__.__name__,
                self.ref_id,
                self.ref_pos,
                self.count_vec
                )

###########################
class CountTab:
    def __init__(self,pos_dict):
        self.pos_dict = pos_dict

    def __getitem__(self,key):
        return self.pos_dict[key]

    def __setitem__(self,key,poscount):
        self.pos_dict.update({key:poscount})

    def __iter__(self):
        return iter(self.pos_dict)

#  construct contig file  #
###########################

class Contig:
    """ Base class for contig """
    def __init__(self, id):
        self.id = id

def initialize_contigs(refpath):
    contigs = {}
    infile = open(refpath,"r" )
    for rec in Bio.SeqIO.parse(infile, 'fasta'):
        contig = Contig(rec.id)
        contig.id = rec.id
        contig.seq = str(rec.seq).upper()
        contig.length = len(contig.seq)
        contigs[contig.id] = contig
    infile.close()
    return contigs

##########################################
#  read bam file and covert to CountTab  #
##########################################

def bam2CountTab(bampath,refpath,baseq,required_coverage):
    bn = os.path.basename(check_path_existence(bampath))
    contigs=initialize_contigs(refpath)
    TableCount = CountTab(pos_dict = OrderedDict())

    pysam.index(bampath)
    with pysam.AlignmentFile(bampath, 'rb') as bamfile:
        for contig_id in sorted(list(contigs.keys())):
            contig = contigs[contig_id]
            num_pos_cov = 0
            counts = bamfile.count_coverage(
                contig.id,
                start=0,
                end=contig.length,
                quality_threshold=baseq,
                read_callback="all")
            covered_bases = 0
            for i in range(0, contig.length):
                check_coverage = counts[0][i] + counts[1][i] + counts[2][i] + counts[3][i]
                if check_coverage > 0:
                    covered_bases += 1
            if covered_bases < contig.length * required_coverage:
                print("The " + str(bampath) + " has less the "+ str(required_coverage*100)+"% genome coverage and it will be removed from analysis")
                return None
            else:
                for i in range(0, contig.length):
                    check_coverage = counts[0][i] + counts[1][i] + counts[2][i] + counts[3][i]
                    if check_coverage > 0:
                        ref_pos = i+1
                        ref_allele = contig.seq[i]
                        count_a = counts[0][i]
                        count_c = counts[1][i]
                        count_g = counts[2][i]
                        count_t = counts[3][i]
                        pos_count = PosCount(
                                ref_id = contig.id,
                                ref_pos = ref_pos,
                                ref_allele = ref_allele,
                                count_vec = np.array([count_a, count_c, count_g, count_t])
                                )
                        TableCount[pos_count.ref_pos] = pos_count
    bamfile.close()
    return (bn, TableCount)

##########################################
#  read tsv files and covert to CountTab  #
##########################################

def rowCount(table_group):
    import pandas as pd
    TableCount = CountTab(pos_dict = OrderedDict())
    for row_index, row in table_group.iterrows():
        pos_count = PosCount(
            ref_id = row["Ref_id"],
            ref_pos = row["Pos"],
            ref_allele = row["Ref"],
            count_vec = std_vec[row["Alt"]])
        TableCount[pos_count.ref_pos] = pos_count
    return TableCount

def table2CountTab(tsv_table):
    import pandas as pd
    #setting up dictionaries
    tsv_dt = OrderedDict()
    #voc_dt = OrderedDict()
    if tsv_table != "":
        #setting up file user input (the user needs to keep the column order)
        tsv = pd.read_csv(tsv_table, sep="\t", header=0, names=["Ref_id", "SNV", "Ref", "Pos", "Alt", "VOC"])
        grouped_tsv = tsv.groupby("VOC")
        for name, tsv_group in grouped_tsv: 
            tsv_tablecount = rowCount(tsv_group)
            tsv_dt[name] = tsv_tablecount
    else:
        None
    return tsv_dt

###############################################
#  construct distance matrix between samples  #
###############################################

def dist2PcoA(dist):
    from skbio.stats.ordination import pcoa
    ordination_result = pcoa(dist)
    return ordination_result


####################################################
#  check the input files and compute the CountTab  #
####################################################

def file2CountTabDict(in_file,refpath,baseq,required_coverage,cpus):
    dt = OrderedDict()
    abf = check_path_existence(in_file)
    bamlist = list()
    with open(abf) as f:
        for each_file in f.read().splitlines():
            if each_file.lower().endswith('.bam'):
                bamlist.append(each_file)
            else:
                sys.exit("File: {} unrecognized file extension! Please check!".format(abspath))
    f.close()

    with Pool(cpus) as p:
        list_dt=p.map(partial(bam2CountTab,refpath=refpath,baseq=baseq,required_coverage=required_coverage),bamlist)
    for item in list_dt:
        if item:
            dt[item[0]]= item[1]
    return dt

def tsvmergetoctd(ctd, tsv_table):
    tsv_ct = table2CountTab(tsv_table)
    dt_merged = {**ctd, **tsv_ct}
    return dt_merged

def dt2dist(dt, cpus):
    from skbio import DistanceMatrix

    key_list = list(dt.keys())
    comb = list(combinations(range(len(dt)), 2))

    inputs_ct = [(dt[key_list[i]],dt[key_list[j]]) for (i,j) in comb]
    with Pool() as p:
        list_dist = p.starmap(CalculateDist, inputs_ct)

    dist = np.zeros(shape=(len(dt),len(dt)))
    idx =0 
    for (i,j) in comb:
        dist[j,i] = list_dist[idx]
        dist[i,j] = dist[j,i]
        idx += 1
    np.fill_diagonal(dist, 0)
    dm = DistanceMatrix(dist,key_list) 
    return dm

def ref2FreqTable(refpath):
    RefTableFreq = CountTab(pos_dict = OrderedDict())
    contigs = initialize_contigs(refpath)
    for contig_id in sorted(list(contigs.keys())):
        contig = contigs[contig_id]
        for pos, base in enumerate(contig.seq):
            if base not in "ACGTacgt":
                base = "N"
            ref_pos_count = PosCount(ref_id = contig.id,
                                 ref_pos = pos + 1,
                                 ref_allele = base,
                                 count_vec = std_vec[base])
            key = ref_pos_count.ref_pos
            RefTableFreq[key] = ref_pos_count
    return(RefTableFreq)

if __name__ == '__main__':

    import logging, time
    from argparse import ArgumentParser
    parser = ArgumentParser(description='Compute pairwise metagenome distance for SARS-CoV-2 samples')
    parser.add_argument('-f', '--file', help='file that lists the path to sorted bam files', required=True, dest='file', metavar='')
    parser.add_argument('-q', '--qual', help='Only reads with mapping quality equal to or greater than this value will be counted (0 by default).', required=False, default=0, dest='baseq', metavar='', type=float)
    parser.add_argument('-c', '--cov', help='Only samples with reads mapped to equal to or greater than this fraction of the genome will be used for PcoA analysis (0.5 by default).', required=False, default=0.5, dest='required_coverage', metavar='', type=float)
    parser.add_argument('-v', '--tsv', help='Option to use a 5 column tab delimited file containing nucleotide mutations from lineages of interest to be added to analysis', default="", required=False, dest='tsv_file', metavar='')
    parser.add_argument('-r', '--ref', help='input reference file', required=False, default="data/NC_045512.2.fasta", dest='refpath', metavar='')
    parser.add_argument('-t', '--threads', help='number of threads used for computing', required=False, default=cpu_count(), dest='cpus', metavar='', type=int)
    parser.add_argument('-o', '--out', help='output folder name', required=True, dest='outpath', metavar='')
    args = parser.parse_args()

    # write log file
    logger = logging.getLogger(__name__)
    logger.setLevel(level=logging.INFO)
    formatter = logging.Formatter('%(asctime)s : %(levelname)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')

    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setLevel(level=logging.INFO)
    logger.addHandler(stream_handler)
    stream_handler.setFormatter(formatter)

    if not os.path.exists(args.outpath):
        os.mkdir(args.outpath)
    else:
        logger.warning('Outpath already exists and will be overwritten!')

    logger.info('Creating count table')
    ctd = file2CountTabDict(args.file, args.refpath,args.baseq,args.required_coverage,args.cpus)
    dt_merged = tsvmergetoctd(ctd, args.tsv_file)

    global RefTableFreq
    RefTableFreq = ref2FreqTable(args.refpath)

    logger.info('Caculating dissimilarity distance')
    dist = dt2dist(dt_merged,args.cpus)
    dist.write("{}/distance.txt".format(args.outpath))

    # logger.info('Obtaining ordination')
    # ordination_result = dist2PcoA(dist)
    # ordination_result.write("{}/ordination_result.txt".format(args.outpath))

    # logger.info('Done!')