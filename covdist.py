#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys, os, subprocess, shutil, csv
import argparse,  platform
import Bio.SeqIO, pysam, numpy as np
from time import time
import pysamstats

from collections import OrderedDict
import numpy as np


std_vec = {
        "A":[1,0,0,0],
        "C":[0,1,0,0],
        "G":[0,0,1,0],
        "T":[0,0,0,1],
        "a":[1,0,0,0],
        "c":[0,1,0,0],
        "g":[0,0,1,0],
        "t":[0,0,0,1]
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

def CalculateDist(count_table1, count_table2):
    sumval = 0
    key_union = set(count_table1).union(count_table2)
    for key in key_union:
        if key in count_table1:
            if key in count_table2:
                sumval += thetayc(count_table1[key].count_vec,count_table2[key].count_vec)
            else:
                sumval += thetayc(count_table1[key].count_vec,std_vec[count_table1[key].ref_allele])
        else:
            sumval += thetayc(count_table2[key].count_vec,std_vec[count_table2[key].ref_allele])
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
    contigs=initialize_contigs(refpath)
    TableCount = CountTab(pos_dict = OrderedDict())

    pysam.index(bampath)
    with pysam.AlignmentFile(bampath, 'rb') as bamfile:
        for contig_id in sorted(list(contigs.keys())):
            contig = contigs[contig_id]
            stat_gen = pysamstats.stat_variation(
                bamfile,
                fafile=refpath,
                chrom=contig.id,
                min_mapq=baseq
                )
            stat_list = list(stat_gen)
            pos_covered = np.count_nonzero([i["reads_all"]>0 for i in stat_list])
            if pos_covered < contig.length * required_coverage:
                print("The " + str(bampath) + " has less the "+ str(required_coverage*100)+"% genome coverage and it will be removed from analysis")
            else:
                for pos in stat_list:
                    if pos["mismatches"] != 0:
                        count_a = pos["A"]
                        count_c = pos["C"]
                        count_g = pos["G"]
                        count_t = pos["T"]
                        pos_count = PosCount(
                                ref_id = contig.id,
                                ref_pos = pos["pos"],
                                ref_allele = pos["ref"],
                                count_vec = np.array([count_a, count_c, count_g, count_t])
                                )
                        TableCount[pos_count.ref_pos] = pos_count

            bamfile.close()
    return TableCount

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

def file2CountTabDict(in_file,refpath,baseq,required_coverage):
    dt = OrderedDict()
    abf = check_path_existence(in_file)
    with open(abf) as f:
        for each_file in f.read().splitlines():
            bn = os.path.basename(check_path_existence(each_file))
#            if bn.lower().endswith(('.fa', '.fasta', '.fna','.fas')):
#                prefix=os.path.splitext(each_file)[0]
#                cmdline = "minimap2 -a {refpath} {each_file}|samtools view -bS|samtools sort > {prefix}.bam".format(refpath=refpath,each_file=each_file,prefix=prefix)
#                os.system(cmdline)
#                dt[bn] = bam2CountTab("{}.bam".format(prefix),refpath)
#            elif bn.lower().endswith('.bam'):
#                dt[bn] = bam2CountTab(each_file,refpath)
            if bn.lower().endswith('.bam'):
                dt[bn] = bam2CountTab(each_file,refpath,baseq,required_coverage)
            else:
                sys.exit("File: {} unrecognized file extension! Please check!".format(abspath))
    f.close()
    return dt

def dt2dist(dt):
    from skbio import DistanceMatrix
    dist = np.zeros(shape=(len(dt),len(dt)))
    key_list = list(dt.keys())
    for i in range(len(dt)):
        for j in range(i):
            dist[j,i] = CalculateDist(dt[key_list[i]],dt[key_list[j]])
            dist[i,j] = CalculateDist(dt[key_list[i]],dt[key_list[j]])
    np.fill_diagonal(dist, 0)
    dm = DistanceMatrix(dist,key_list) 
    return dm


if __name__ == '__main__':

    import logging, time
    from argparse import ArgumentParser
    parser = ArgumentParser(description='Compute pairwise metagenome distance for SARS-CoV-2 samples')
    parser.add_argument('-f', '--file',help='file that lists the path to sorted bam files', required=True,dest='file',metavar='')
    parser.add_argument('-m', '--meta',help='file containing metadata info', required=True,dest='meta',metavar='')
    parser.add_argument('-q', '--qual',help='Only reads with mapping quality equal to or greater than this value will be counted (0 by default).', required=False,default=0,dest='baseq',metavar='')
    parser.add_argument('-c', '--cov',help='Only samples with reads mapped to equal to or greater than this fraction of the genome will be used for PcoA analysis (0.5 by default).', required=False,default=0.5,dest='required_coverage',metavar='')
    parser.add_argument('-r', '--ref', help='input reference file', required=False,default="data/NC_045512.2.fasta",dest='refpath',metavar='')
    parser.add_argument('-o', '--out', help='output folder name', required=True,dest='outpath',metavar='')
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
    ctd = file2CountTabDict(args.file, args.refpath,args.baseq,args.required_coverage)
    logger.info('Caculating dissimilarity distance')
    dist = dt2dist(ctd)
    dist.write("{}/distance.txt".format(args.outpath))

    logger.info('Obtaining ordination')
    ordination_result = dist2PcoA(dist)
    ordination_result.write("{}/ordination_result.txt".format(args.outpath))

    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns

    meta = pd.read_csv(args.meta,sep="\t",index_col=0,header=None,names=["info"])
    meta['info']=meta['info'].astype("category")
    df = pd.concat([ordination_result.samples,meta],axis=1)

    logger.info('Generate PcoA plots')
    sns.lmplot( x="PC1", y="PC2", data=df, fit_reg=False, hue='info', legend=False, palette="Set1")
    plt.legend(loc='lower right')
    plt.gca().set_aspect("equal",adjustable="box")
    plt.savefig("{}/PcoA_plot.svg".format(args.outpath))
    plt.savefig("{}/PcoA_plot.png".format(args.outpath))

    plt2 = ordination_result.plot( df=meta, column='info',cmap='Set1', s=50)
    plt2.savefig("{}/PcoA_plot_3D.svg".format(args.outpath))
    plt2.savefig("{}/PcoA_plot_3D.png".format(args.outpath))
