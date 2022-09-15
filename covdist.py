#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys, os
import logging, vcf, click
import Bio.SeqIO, numpy as np
from multiprocessing import Pool, cpu_count
from functools import partial
from itertools import combinations
from collections import OrderedDict
import pandas as pd

std_vec = {
        "A":[1.0,0.0,0.0,0.0],
        "C":[0.0,1.0,0.0,0.0],
        "G":[0.0,0.0,1.0,0.0],
        "T":[0.0,0.0,0.0,1.0],
        "a":[1.0,0.0,0.0,0.0],
        "c":[0.0,1.0,0.0,0.0],
        "g":[0.0,0.0,1.0,0.0],
        "t":[0.0,0.0,0.0,1.0],
        "N":[0.0,0.0,0.0,0.0]
        }

####################
#  define classes  #
####################

class PosFreq:
    def __init__(self,ref_id,ref_pos,ref_allele,freq_vec):
        self.ref_id = ref_id
        self.ref_pos = ref_pos
        self.ref_allele = ref_allele
        self.freq_vec = freq_vec

    def __repr__(self):
        # !r means use __repr__ insteadd of __str__
        return "{0}(id={1!r}, pos={2!r}, freq_vec={3!r})".format(
                self.__class__.__name__,
                self.ref_id,
                self.ref_pos,
                self.freq_vec
                )
                
class FreqTab:
    def __init__(self,pos_dict):
        self.pos_dict = pos_dict
        self.filtered_pos = []

    def __getitem__(self,key):
        return self.pos_dict[key]

    def __setitem__(self,key,posfreq):
        self.pos_dict.update({key:posfreq})

    def __iter__(self):
        return iter(self.pos_dict)
    
    def __repr__(self):
        string = ""
        for key in self.pos_dict:
            string += str(key) + '_' + repr(self.pos_dict[key]) + "\n"
        return(string)

    def filter_pos(self, pos_depth_dict, depth_cutoff=0):
        # exclude low depth positions
        for key in pos_depth_dict:
            if pos_depth_dict[key] <= depth_cutoff:
                # del self.pos_dict[key]
                self.filtered_pos.append(key)

class Contig:
    """ Base class for contig """
    def __init__(self, id):
        self.id = id

#########################
#  essential functions  #
#########################

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

def check_path_existence(file_path):
    abspath = os.path.abspath(file_path)
    if os.path.exists(abspath):
        return(abspath)
    else:
        return(None)

def calculate_coverage(pos_depth_dict, contig, cov_depth_cutoff=10):
    base_cov = 0
    for key in pos_depth_dict:
        ref_id = key.split(':')[0]
        if ref_id == contig.id:
            if pos_depth_dict[key] >= cov_depth_cutoff:
                base_cov += 1
    return base_cov/contig.length

def thetayc(vec1, vec2):
    # calculate Yue & Clayton measure of dissimilarity
    if np.sum(vec1) == 0 or np.sum(vec2) ==0:
        return 0
    else:
        dotsum = np.dot(vec1, vec2)
        dist = np.sum(np.square(vec1),axis=0) + np.sum(np.square(vec2),axis=0) 
        if dist == dotsum:
            return 0
        else:
            return 1-dotsum/(dist-dotsum)

# def round_vec(vec):
#     return [round(i,3) for i in vec]

def CalculateDist(freq_table1, freq_table2):
    ### debug ###
    # with open("test1", "w") as g:
    #     for key in RefTableFreq:
    #         if key in freq_table1:
    #             vec1 = round_vec(freq_table1[key].freq_vec)
    #         else:
    #             vec1 = round_vec(RefTableFreq[key].freq_vec)
    #         if key in freq_table2:
    #             vec2 = round_vec(freq_table2[key].freq_vec)
    #         else:
    #             vec2 = round_vec(RefTableFreq[key].freq_vec)
    #         val = thetayc(vec1, vec2)
    #         v1 = "\t".join(list(map(str, vec1)))
    #         v2 = "\t".join(list(map(str, vec2)))
    #         g.write(str(key.split(":")[1])+"\t"+v1+"\t"+v2+"\t"+str(val)+"\n")

    sumval = 0
    key_inter = set(freq_table1).intersection(set(freq_table2))
    for key in key_inter:
        if key not in freq_table1.filtered_pos + freq_table2.filtered_pos:
            sumval += thetayc(freq_table1[key].freq_vec, freq_table2[key].freq_vec)
    for key in (set(freq_table1) - key_inter):
        if key not in freq_table1.filtered_pos + freq_table2.filtered_pos:
            sumval += thetayc(freq_table1[key].freq_vec, RefTableFreq[key].freq_vec)
    for key in (set(freq_table2) - key_inter):
        if key not in freq_table1.filtered_pos + freq_table2.filtered_pos:
            sumval += thetayc(freq_table2[key].freq_vec, RefTableFreq[key].freq_vec)
    return sumval

#############################################
#  read depth file and covert to dictionary #
#############################################

def depth2Dict(depth_file):
    pos_depth_dict = OrderedDict()
    df = pd.read_csv(depth_file, sep="\t", header=0, names=["Ref_id", "Pos", "Ref", "Depth"])
    for _, row in df.iterrows():
        pos_depth_dict[row['Ref_id']+":"+str(row['Pos'])] = int(row['Depth'])
    return pos_depth_dict

##########################################
#  read ref file and covert to FreqTab   #
##########################################

def ref2FreqTable(refpath):
    RefTableFreq = FreqTab(pos_dict = OrderedDict())
    contigs = initialize_contigs(refpath)
    for contig_id in sorted(list(contigs.keys())):
        contig = contigs[contig_id]
        for pos, base in enumerate(contig.seq):
            if base not in "ACGTacgt":
                base = "N"
            ref_pos_freq = PosFreq(ref_id = contig.id,
                                 ref_pos = pos + 1,
                                 ref_allele = base,
                                 freq_vec = std_vec[base])
            key = ref_pos_freq.ref_id + ":" + str(ref_pos_freq.ref_pos)
            RefTableFreq[key] = ref_pos_freq
    return(RefTableFreq)

##########################################
#  read vcf file and covert to FreqTab   #
##########################################

def vcf2FreqTable(prefix, refpath, cov_depth_cutoff, required_coverage):
    # check vcf file
    vcf_abf = check_path_existence(prefix+'.vcf')
    if not vcf_abf:
        sys.exit("File: {} doesn't exist. Please check!".format(vcf_abf))

    # check depth file
    depth_abf = check_path_existence(prefix+'.depth')
    if not depth_abf:
        logger.warning('{} is not found. Default to set all position depths are valid.'.format(depth_abf))

    prefix = os.path.basename(vcf_abf).replace(".vcf", "")
    TableFreq = FreqTab(pos_dict = OrderedDict())
    
    contigs = initialize_contigs(refpath)

    for contig_id in sorted(list(contigs.keys())):
        contig = contigs[contig_id]
        vcf_reader = vcf.Reader(open(vcf_abf, 'r'))
        for record in vcf_reader:
            if len(record.FILTER) >= 1 and record.FILTER[0] == 'AmpliconRemoval':
                continue
            if record.is_snp and record.CHROM == contig.id:
                freq_vector = [0.0, 0.0, 0.0, 0.0]
                for index, alt in enumerate(record.ALT):
                    if str(alt) in 'ACGT':
                        if 'AF' in record.INFO:
                            if type(record.INFO['AF']) != float:
                                freq_vector['ACGT'.index(str(alt))] = record.INFO['AF'][index]
                            else:
                                freq_vector['ACGT'.index(str(alt))] = record.INFO['AF']
                        else:
                            freq_vector['ACGT'.index(str(alt))] = 1
                
                # get reference base frequency
                if 1 - sum(freq_vector) > 1e-5:
                    freq_vector['ACGT'.index(str(record.REF))] = 1 - sum(freq_vector)

                pos_freq = PosFreq(
                    ref_id = record.CHROM,
                    ref_pos = record.POS,
                    ref_allele = record.REF,
                    freq_vec = freq_vector)
            key = pos_freq.ref_id + ":" + str(pos_freq.ref_pos)

            TableFreq[key] = pos_freq
        
        if depth_abf:
            # if depth path is provided, filter positions be depth
            pos_depth_dict = depth2Dict(depth_abf)
            TableFreq.filter_pos(pos_depth_dict, depth_cutoff=0)

            coverage = calculate_coverage(pos_depth_dict, contig, cov_depth_cutoff)
            if coverage < required_coverage:
                logger.warning('{0} coverage ({1}%) is less than {2}% genome coverage and it will be removed from analysis!'.format(vcf_abf, coverage*100, required_coverage*100))

    return (prefix, TableFreq)

###############################################
#  parse vcf files to a dict of freq tables   #
###############################################

def vcf2FreqTabDict(prefix_list, refpath, cov_depth_cutoff, required_coverage, cpus):
    dt = OrderedDict()
    prefix_abf = check_path_existence(prefix_list)
    if prefix_abf:
        prefixlist = [line.strip() for line in open(prefix_abf)]
    else:
        sys.exit("File: {} doesn't exist. Please check!".format(prefix_abf))
    
    with Pool(cpus) as p:
        list_dt = p.map(partial(vcf2FreqTable, refpath=refpath, cov_depth_cutoff=cov_depth_cutoff, required_coverage=required_coverage), prefixlist)

    for item in list_dt:
        if item:
            dt[item[0]]= item[1]
    return dt

###############################################
#  construct distance matrix between samples  #
###############################################

def dt2dist(dt, cpus):
    from skbio import DistanceMatrix

    key_list = list(dt.keys())
    comb = list(combinations(range(len(dt)), 2))

    inputs_ct = [(dt[key_list[i]], dt[key_list[j]]) for (i,j) in comb]

    with Pool(cpus) as p:
        list_dist = p.starmap(CalculateDist, inputs_ct)

    dist = np.zeros(shape=(len(dt),len(dt)))
    idx = 0
    for (i,j) in comb:
        dist[j,i] = list_dist[idx]
        dist[i,j] = dist[j,i]
        idx += 1

    np.fill_diagonal(dist, 0)
    dm = DistanceMatrix(dist,key_list) 
    return(dm)

############################################
#  convert distance matrix to PcoA result  #
############################################

def dist2PcoA(dist):
    from skbio.stats.ordination import pcoa
    ordination_result = pcoa(dist)
    return ordination_result

####################################
###  main function for cov-dist  ###
####################################

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
@click.group(chain=True, context_settings=CONTEXT_SETTINGS)
def cli():
    global logger
    logger = logging.getLogger(__name__)
    logger.setLevel(level=logging.INFO)
    formatter = logging.Formatter('%(asctime)s : %(levelname)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')

    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setLevel(level=logging.INFO)
    logger.addHandler(stream_handler)
    stream_handler.setFormatter(formatter)

@cli.command()
@click.option('-p', '--prefix_list', required=True, type=click.Path(exists=True), help='A file that lists the prefix for vcf and corresponidng depth files.')
@click.option('-r', '--ref', 'refpath', required=False, default=os.path.join(os.path.abspath(os.path.dirname( __file__ )), "data", "NC_045512.2.fasta"), type=click.Path(exists=True), help='Input reference file. [default: data/NC_045512.2.fasta]')
@click.option('-c', '--cov', 'required_coverage', required=False, default=0.5, type=float, help='Only samples with reads mapped to equal to or greater than this fraction of the genome will be used for PcoA analysis [default: 0.5].')
@click.option('-d', '--depth', 'cov_depth_cutoff', required=False, default=10, type=int, help='Depth cutoff to include positions when calculating coverage [default: 10].')
@click.option('-t', '--threads', 'cpus', required=False, default=cpu_count(), type=int, help='Number of threads used for computing [default: all available cpus].')
# @click.option('-v', '--voc_list', 'tsv_file', required=False, default="", help='Option to use a 5 column tab delimited file containing nucleotide mutations from lineages of interest to be added to analysis')
@click.option('-o', '--out', 'outpath', required=True, help='Output folder name')
def dist(prefix_list, refpath, cov_depth_cutoff, required_coverage, cpus, outpath):
    """Compute pairwise metagenome distance for SARS-CoV-2 samples."""
    if not os.path.exists(outpath):
        os.mkdir(outpath)
    else:
        logger.warning('Outpath already exists and will be overwritten!')

    logger.info('Creating frequency table')
    freq_tab_dict = vcf2FreqTabDict(prefix_list, refpath, cov_depth_cutoff, required_coverage, cpus)
    global RefTableFreq
    RefTableFreq = ref2FreqTable(refpath)

    logger.info('Caculating dissimilarity distance')
    dist = dt2dist(freq_tab_dict, cpus)
    dist.write("{}/distance.txt".format(outpath))
    
    logger.info('Calculating distances is done!')

@cli.command()
@click.option('-d', '--distance_file', type=click.Path(exists=True), required=True, help='Distance file obtained from cov-dist "dist" command.')
@click.option('-m', '--meta', type=click.Path(exists=True), required=True, help='Meta data for samples.')
@click.option('-c', '--column', type=str, required=True, help='The column name in the meta data to color the samples.')
@click.option('-o', '--out', 'outpath', required=True, help='Output folder name')
def plot(distance_file, meta, column, outpath):
    """Perform ordination analysis and visualize results."""
    import pandas as pd
    import plotly.express as px

    if not os.path.exists(outpath):
        os.mkdir(outpath)

    logger.info('Obtaining ordination')
    dist = pd.read_csv(distance_file, sep="\t", header=0, index_col=0)
    ordination_result = dist2PcoA(dist)
    ordination_result.write("{}/ordination_result.txt".format(outpath))

    logger.info('Performing ordination analysis is done!')

    logger.info('Plotting')
    metadata = pd.read_csv(meta, sep="\t", header=0)
    if "sample" not in metadata.columns:
        logger.error('Column "sample" is not in meta data. Please assign a column to match the columns in distance matrix.')
    if column not in metadata.columns:
        logger.error('Column "{}" is not in meta data. Please check!'.format(column))

    sample_names = list(dist.index)
    if len(sample_names) >= 3:
        dt = ordination_result.samples.loc[:,['PC1', 'PC2', 'PC3']]
        dt['sample'] = sample_names
        plot_dt = dt.join(metadata.set_index('sample'), on='sample')
        fig = px.scatter_3d(plot_dt, x='PC1', y='PC2', z='PC3', color=column)
        fig.write_html("{}/pcoa_plot.html".format(outpath))
    elif len(sample_names) == 2:
        dt = ordination_result.samples.loc[:,['PC1', 'PC2']]
        dt['sample'] = sample_names
        plot_dt = dt.join(metadata.set_index('sample'), on='sample')
        fig = px.scatter(plot_dt, x='PC1', y='PC2', color=column)
        fig.write_html("{}/pcoa_plot.html".format(outpath))
    else:
        logger.warning('Only one sample is found. No plot is generated.')
    
    logger.info('Plotting is done!')

if __name__ == '__main__':
    cli()