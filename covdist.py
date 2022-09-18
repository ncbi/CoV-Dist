#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys, os, datetime, logging, vcf, click
import Bio.SeqIO
import numpy as np
from multiprocessing import Pool, cpu_count
from functools import partial
from itertools import combinations
from collections import OrderedDict
import pandas as pd
import plotly.express as px
import plotly.graph_objs as go

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

def CalculateDist(freq_table1, freq_table2):
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

def vcf2FreqTable(prefix, refpath, cov_depth_cutoff, required_coverage, mode="normal"):
    
    prefix_base = os.path.basename(prefix)

    # check vcf file
    vcf_abf = check_path_existence(prefix+'.vcf')
    if not vcf_abf:
        raise FileNotFoundError("Vcf file: {} doesn't exist! Terminate program.".format(prefix+'.vcf'))

    # check depth file
    depth_abf = check_path_existence(prefix+'.depth')

    if not depth_abf and mode == "normal":
        logger.warning('Depth file: {} is not found. Default to set all position depths are valid.'.format(prefix+'.depth'))
    
    TableFreq = FreqTab(pos_dict = OrderedDict())
    contigs = initialize_contigs(refpath)
    for contig_id in sorted(list(contigs.keys())):
        contig = contigs[contig_id]

        if depth_abf:
            # if depth path is provided, filter positions be depth
            pos_depth_dict = depth2Dict(depth_abf)
            TableFreq.filter_pos(pos_depth_dict, depth_cutoff=0) # this can be modified
            coverage = calculate_coverage(pos_depth_dict, contig, cov_depth_cutoff)
            if coverage < required_coverage:
                logger.warning('{0} coverage ({1}%) is less than {2}% genome coverage and it will be removed from analysis!'.format(vcf_abf, coverage*100, required_coverage*100))
                return None

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

    return (prefix_base, TableFreq)

###############################################
#  parse vcf files to a dict of freq tables   #
###############################################

def vcf2FreqTabDict(prefixlist, refpath, cov_depth_cutoff, required_coverage, cpus, mode="normal"):
    dt = OrderedDict()
    
    with Pool(cpus) as p:
        try:
            list_dt = p.map(partial(vcf2FreqTable, refpath=refpath, cov_depth_cutoff=cov_depth_cutoff, required_coverage=required_coverage, mode=mode), prefixlist)
        except FileNotFoundError as e:
            logger.error(e)
            sys.exit()

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

    with Pool(cpus) as p1:
        list_dist = p1.starmap(CalculateDist, inputs_ct)

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
@click.group(context_settings=CONTEXT_SETTINGS)
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
@click.option('-v', '--voc', 'voc_dir', required=False, default=None, help='voc folder to store vcf files (default is None).')
@click.option('-o', '--out', 'outpath', required=True, help='Output folder name')
def dist(prefix_list, refpath, cov_depth_cutoff, required_coverage, cpus, voc_dir, outpath):
    """Compute pairwise metagenome distance for SARS-CoV-2 samples."""
    if not os.path.exists(outpath):
        os.mkdir(outpath)
    else:
        logger.warning('Outpath already exists and will be overwritten!')

    logger.info('Creating frequency table')

    prefixlist = [line.strip() for line in open(prefix_list)]

    freq_tab_dict = vcf2FreqTabDict(prefixlist, refpath, cov_depth_cutoff, required_coverage, cpus, mode="normal")
    global RefTableFreq
    RefTableFreq = ref2FreqTable(refpath)

    if not voc_dir is None:
        logger.info('Creating voc frequency table')
        voc_prefixlist = []
        for file in os.listdir(voc_dir):
            voc_prefix = os.path.join(voc_dir, file).replace(".vcf", "")
            voc_prefixlist.append(voc_prefix)
        voc_freq_tab_dict = vcf2FreqTabDict(voc_prefixlist, refpath, cov_depth_cutoff, required_coverage, cpus, mode="voc")
        final_freq_tab_dict = {**freq_tab_dict, **voc_freq_tab_dict}
    else:
        final_freq_tab_dict = freq_tab_dict

    logger.info('Calculating dissimilarity distance')
    dist = dt2dist(final_freq_tab_dict, cpus)
    dist.write("{}/distance.txt".format(outpath))
    
    logger.info('Calculating distances is done!')

###########################
# plot essential funtions #
###########################
my_symbols = ['diamond', 'x', 'circle', 'square', 'cross', 'diamond-open', 'circle-open', 'square-open']

def validate_date_format(date_text):
    try:
        datetime.datetime.strptime(date_text, '%Y-%m-%d')
    except ValueError:
        raise ValueError('Incorrect data format for "{}", should be YYYY-MM-DD!'.format(date_text))

def generate_hovertemplate(df_columns, method="PCoA"):
    if method == "PCoA":
        pc_name = "PC"
    elif method == "MDS":
        pc_name = "Dimension"
    elif method == "TSNE":
        pc_name = "t-SNE"

    txt_list = ["sample: %{customdata[0]}",
            "collection_date: %{customdata[1]}",
            "collection_site: %{customdata[2]}",
            pc_name+"1: %{customdata[3]}",
            pc_name+"2: %{customdata[4]}",
            pc_name+"3: %{customdata[5]}"
            ]
    for index, col in enumerate(df_columns):
        if col not in ["sample", "collection_date", "collection_site", "PC1", "PC2", "PC3"]:
            txt_list += [str(col)+": %{customdata["+str(index)+"]}"]
    
    return("<br>".join(txt_list))

def check_meta_column(column, target_columns):
    if column !="default" and column not in target_columns:
        logger.error('Column "{}" is not in meta data. Please check!'.format(column))
        sys.exit()

def other_col_plot(df, column, analysis, intersect_len, outpath, axis_names, hovertmp, voc=False, mode="3d"):
    
    fig = go.Figure()
    
    if mode == "3d":
        if voc == False:
            plot3d = px.scatter_3d(df, x='PC1', y='PC2', z='PC3', color=column, color_continuous_scale="Viridis", hover_data=list(df.columns))
            fig.add_traces(list(plot3d.select_traces()))
        else:
            df1 = df[df['voc']=='no']
            plot3d = px.scatter_3d(df1, x='PC1', y='PC2', z='PC3', color=column, color_continuous_scale="Viridis", hover_data=list(df1.columns))
            fig.add_traces(list(plot3d.select_traces()))

            df2 = df[df['voc']=='yes']
            unique_type = df2['type'].unique()
            for index, typ in enumerate(unique_type):
                df_typ = df2[df2['type'] == typ]
                legendgroup=None
                legendgrouptitle_text=None
                if index==0:
                    legendgroup="clade"
                    legendgrouptitle_text="<br>clade"
                trace = go.Scatter3d(x=df_typ['PC1'], y=df_typ['PC2'], z=df_typ['PC3'], 
                                mode='markers',
                                marker_size=6,
                                marker_color="#656565",
                                marker_showscale=False,
                                marker_symbol=[my_symbols[index%len(my_symbols)]]*df_typ.shape[0],
                                name=typ, 
                                hoverinfo='skip', # no hover text
                                showlegend=True,
                                customdata=df2,
                                hovertemplate=hovertmp,
                                legendgroup=legendgroup,
                                legendgrouptitle_text=legendgrouptitle_text
                )
                fig.add_trace(trace)

        fig.update_layout(title="CoV-Dist {0} {1} with {2} samples".format(mode, analysis, intersect_len), 
                    scene = dict(
                    xaxis_title=axis_names[0],
                    yaxis_title=axis_names[1],
                    zaxis_title=axis_names[2]))

    elif mode == "2d":
        if voc == False:
            plot2d = px.scatter(df, x='PC1', y='PC2', color=column, color_continuous_scale="Viridis", hover_data=list(df.columns))
            fig.add_traces(list(plot2d.select_traces()))
        else:
            df1 = df[df['voc']=='no']
            plot2d = px.scatter(df1, x='PC1', y='PC2', color=column, color_continuous_scale="Viridis", hover_data=list(df1.columns))
            fig.add_traces(list(plot2d.select_traces()))

            df2 = df[df['voc']=='yes']
            unique_type = df2['type'].unique()
            for index, typ in enumerate(unique_type):
                df_typ = df2[df2['type'] == typ]
                legendgroup=None
                legendgrouptitle_text=None
                if index==0:
                    legendgroup="clade"
                    legendgrouptitle_text="<br>clade"
                trace = go.Scatter(x=df_typ['PC1'], y=df_typ['PC2'], 
                                mode='markers',
                                marker_size=6,
                                marker_color="#656565",
                                marker_showscale=False,
                                marker_symbol=[my_symbols[index%len(my_symbols)]]*df_typ.shape[0],
                                name=typ, 
                                customdata=df2,
                                hovertemplate=hovertmp,
                                showlegend=True,
                                legendgroup=legendgroup,
                                legendgrouptitle_text=legendgrouptitle_text
                )
                fig.add_trace(trace)

        fig.update_layout(title="CoV-Dist {0} {1} with {2} samples".format(mode, analysis, intersect_len), 
                    xaxis_title=axis_names[0],
                    yaxis_title=axis_names[1])
    
    fig.update_layout(legend_title_text=column)
    
    fig.update_layout(legend=dict(
        orientation="v",
        yanchor="top",
        y=1,
        xanchor="right",
        x=-0.1,
        tracegroupgap=3
    ))

    fig.write_html("{0}/{1}_plot_{2}.html".format(outpath, analysis.lower(), mode))

def default_plot(df, analysis, intersect_len, outpath, axis_names, hovertmp, voc=False, mode="3d"):
    
    if voc == False:
        df_non_voc = df[:]
    else:
        df_non_voc = df[df['voc']=='no'][:]

    df_non_voc['collection_date']=pd.to_datetime(df_non_voc['collection_date']).dt.date
    df_non_voc = df_non_voc.sort_values(by='collection_date')
    unique_date = df_non_voc['collection_date'].unique()
    date_to_val = df_non_voc['collection_date'].map(pd.Series(data=np.arange(len(unique_date)), index=unique_date).to_dict())
    df_non_voc['date2val'] = date_to_val

    step = max(1, int((max(date_to_val)-min(date_to_val))/10))
    tickvals = [k for k in range(min(date_to_val), max(date_to_val), step)]
    dlist = list(date_to_val)
    index_tickvals = [dlist.index(tv) for tv in tickvals]
    ticktext = [df_non_voc['collection_date'][id].strftime("%Y-%m-%d") for id in index_tickvals]
    unique_loc = df_non_voc['collection_site'].unique()
    df = df.join(df_non_voc[['sample', 'date2val']].set_index('sample'), on='sample', how='left')

    fig = go.Figure()
    if mode=="3d":
        if voc == False:
            for loc in unique_loc:
                df_loc = df[df['collection_site'] == loc]
                trace = go.Scatter3d(x=df_loc['PC1'], y=df_loc['PC2'], z=df_loc['PC3'], 
                                mode='markers',
                                marker_size=6,
                                marker_color=df_loc['date2val'],
                                marker_colorscale='Plasma',
                                marker_showscale=False,
                                name=loc, 
                                showlegend=True)
                fig.add_trace(trace)
        else:
            df1 = df[df['voc']=='no']
            for index, loc in enumerate(unique_loc):
                df_loc = df1[df1['collection_site'] == loc]
                trace = go.Scatter3d(x=df_loc['PC1'], y=df_loc['PC2'], z=df_loc['PC3'], 
                                mode='markers',
                                marker_size=6,
                                marker_color=df_loc['date2val'],
                                marker_colorscale='Plasma',
                                marker_showscale=False,
                                name=loc,
                                showlegend=True)
                fig.add_trace(trace)

            df2 = df[df['voc']=='yes']
            unique_type = df2['type'].unique()
            for index, typ in enumerate(unique_type):
                df_typ = df2[df2['type'] == typ]
                legendgroup=None
                legendgrouptitle_text=None
                if index==0:
                    legendgroup="clade"
                    legendgrouptitle_text="<br>clade"
                trace = go.Scatter3d(x=df_typ['PC1'], y=df_typ['PC2'], z=df_typ['PC3'], 
                                mode='markers',
                                marker_size=6,
                                marker_color="#656565",
                                marker_showscale=False,
                                marker_symbol=[my_symbols[index % len(my_symbols)]]*df_typ.shape[0],
                                name=typ, 
                                hoverinfo='skip', # no hover text
                                showlegend=True,
                                legendgroup=legendgroup,
                                legendgrouptitle_text=legendgrouptitle_text)

                fig.add_trace(trace)
        
        df = df.drop('date2val', axis=1)

        fig.add_trace(go.Scatter3d(x=df['PC1'], y=df['PC2'], z=df['PC3'], 
                        mode='markers',
                        marker_color=date_to_val,
                        marker_colorscale='Plasma',
                        marker_showscale=True,
                        marker_size=6,
                        name="",
                        opacity=0,
                        showlegend=False,
                        marker_colorbar=dict(tickvals=tickvals, ticktext=ticktext, title_text='collection_date'),
                        customdata=df,
                        hoverlabel = dict(namelength=0),
                        hovertemplate=hovertmp))

        fig.update_layout(title="CoV-Dist {0} {1} with {2} samples".format(mode, analysis, intersect_len), 
            scene = dict(
            xaxis_title=axis_names[0],
            yaxis_title=axis_names[1],
            zaxis_title=axis_names[2]))
    
    elif mode=="2d":
        if voc == False:
            for loc in unique_loc:
                df_loc = df[df['collection_site'] == loc]
                # add collection_site trace
                trace = go.Scatter(x=df_loc['PC1'], y=df_loc['PC2'], 
                                mode='markers',
                                marker_size=6,
                                marker_color=df_loc['date2val'],
                                marker_colorscale='Plasma',
                                marker_showscale=False,
                                name=loc, 
                                showlegend=True)
                fig.add_trace(trace)
        else:
            df1 = df[df['voc']=='no']
            for index, loc in enumerate(unique_loc):
                df_loc = df1[df1['collection_site'] == loc]
                trace = go.Scatter(x=df_loc['PC1'], y=df_loc['PC2'], 
                                mode='markers',
                                marker_size=6,
                                marker_color=df_loc['date2val'],
                                marker_colorscale='Plasma',
                                marker_showscale=False,
                                name=loc,
                                showlegend=True)
                fig.add_trace(trace)

            df2 = df[df['voc']=='yes']
            unique_type = df2['type'].unique()
            for index, typ in enumerate(unique_type):
                df_typ = df2[df2['type'] == typ]
                legendgroup=None
                legendgrouptitle_text=None
                if index==0:
                    legendgroup="clade"
                    legendgrouptitle_text="<br>clade"
                trace = go.Scatter(x=df_typ['PC1'], y=df_typ['PC2'], 
                                mode='markers',
                                marker_size=6,
                                marker_color="#656565",
                                marker_showscale=False,
                                marker_symbol=[my_symbols[index % len(my_symbols)]]*df_typ.shape[0],
                                name=typ, 
                                hoverinfo='skip', # no hover text
                                showlegend=True,
                                legendgroup=legendgroup,
                                legendgrouptitle_text=legendgrouptitle_text)

                fig.add_trace(trace)

        df = df.drop('date2val', axis=1)
        fig.add_trace(go.Scatter(x=df['PC1'], y=df['PC2'], 
                        mode='markers',
                        marker_color=date_to_val,
                        marker_colorscale='Plasma',
                        marker_showscale=True,
                        marker_size=6,
                        name="",
                        opacity=0,
                        marker_colorbar=dict(tickvals=tickvals, ticktext=ticktext, title_text='collection_date'),
                        customdata=df,
                        hoverlabel = dict(namelength=0),
                        hovertemplate=hovertmp))
        
        fig.update_layout(title="CoV-Dist {0} {1} with {2} samples".format(mode, analysis, intersect_len), 
                    xaxis_title=axis_names[0],
                    yaxis_title=axis_names[1])

    # both for 3d and 2d
    fig.update_layout(legend_title_text="collection_site")

    fig.update_layout(legend=dict(
        orientation="v",
        yanchor="top",
        y=1,
        xanchor="right",
        x=-0.1,
        tracegroupgap=3
    ))

    fig.write_html("{0}/{1}_plot_{2}.html".format(outpath, analysis.lower(), mode))

@cli.command()
@click.option('-d', '--distance_file', type=click.Path(exists=True), required=True, help='Distance file obtained from cov-dist "dist" command.')
@click.option('-m', '--meta', type=click.Path(exists=True), required=True, help='Meta data for samples. Must have at least three columns - "sample", "collection_date" and "collection_site".')
@click.option('-c', '--column', type=str, default="default", required=False, help='The column name in the meta data to color the samples. Default collection_date and collection_site are used to plot figures.')
@click.option('-a', '--analysis', type=click.Choice(['PCoA', 'MDS', 'TSNE'], case_sensitive=False), required=False, default="PCoA", help='Method to show samples. Can be choosen from PCoA (default), MDS, TSNE.')
@click.option('-v', '--voc_meta', type=click.Path(exists=True), required=False, default=None, help='Metadata for VOC.')
@click.option('-o', '--out', 'outpath', required=True, help='Output folder name')
def plot(distance_file, meta, column, analysis, voc_meta, outpath):
    """Perform ordination analysis and visualize results."""

    import warnings
    warnings.simplefilter(action='ignore', category=FutureWarning)
    warnings.simplefilter(action='ignore', category=RuntimeWarning)

    if not os.path.exists(outpath):
        os.mkdir(outpath)

    logger.info('Performing {} ordination analysis'.format(analysis))
    dist = pd.read_csv(distance_file, sep="\t", header=0, index_col=0)
    sample_names = list(dist.index.astype('string'))
    if len(sample_names) < 3:
        logger.error('Less than 3 samples are found. No plot is generated.')
        sys.exit()

    # read metadata
    metadata = pd.read_csv(meta, sep="\t", header=0)
    # check sample column
    check_meta_column("sample", metadata.columns)

    # check collection_date column
    check_meta_column("collection_date", metadata.columns)
    for date_str in metadata['collection_date']:
        try:
            validate_date_format(date_str)
        except ValueError as e:
            logger.error(e)
            sys.exit()
    # check collection_site column
    check_meta_column("collection_site", metadata.columns)
    # check user custom column
    check_meta_column(column, metadata.columns)

    if voc_meta is not None:
        metadt_voc = pd.read_csv(voc_meta, sep="\t", header=0)
        check_meta_column("sample", metadt_voc.columns)
        check_meta_column("type", metadt_voc.columns)
        metadata['voc']='no'
        metadt_voc['voc']='yes'
        merged_metadata = pd.concat([metadata, metadt_voc], ignore_index=True)
        voc = True
    else:
        merged_metadata = metadata
        voc = False

    merged_metadata['sample'] = merged_metadata['sample'].astype('string')
    merged_metadata['collection_site'] = merged_metadata['collection_site'].astype('string')
    merged_metadata['collection_date'] = pd.to_datetime(merged_metadata['collection_date']).dt.date

    # check if distance samples are in metadata
    diff_samples = set(sample_names).difference(set(list(merged_metadata['sample'])))
    diff_len = len(diff_samples)
    if diff_len == 1:
       logger.warning('One sample ({}) is not found in metadata.'.format(diff_samples[0]))
    elif 1 < diff_len <= 10:
        logger.warning('{0} samples ({1}) are not found in metadata.'.format(diff_len, ", ".join(list(diff_samples))))
    elif diff_len > 10:
        logger.warning('{0} samples ({1}) ... are not found in metadata.'.format(diff_len, ", ".join(list(diff_samples)[:10])))
    
    intersect_samples = set(sample_names).intersection(set(list(merged_metadata['sample'])))
    intersect_len = len(intersect_samples)
    if intersect_len == 0:
        logger.error('No samples are not found in metadata.')
        sys.exit()

    ########
    # PcoA #
    ########
    if analysis == "PCoA":
        ordination_result = dist2PcoA(dist)
        dt = ordination_result.samples.loc[:,['PC1', 'PC2', 'PC3']]
        dt['sample'] = sample_names
        axis_names = ["PC" + str(index+1) + " ("+ '{:.2f}'.format(prop*100) +"%)" for index, prop in enumerate(ordination_result.proportion_explained[:3])]

    #######
    # MDS #
    #######
    elif analysis == "MDS":
        from sklearn.manifold import MDS
        mds = MDS(n_components=3, dissimilarity="precomputed")
        mds_fit = mds.fit(dist)
        dt = pd.DataFrame(mds_fit.embedding_, columns=['PC1', 'PC2', 'PC3'])
        dt['sample'] = sample_names
        axis_names = ["Dimension 1", "Dimension 2", "Dimension 3"]

    #########
    # t-SNE #
    #########
    elif analysis == "TSNE":
        from sklearn.manifold import TSNE
        tsne = TSNE(n_components=3, metric="precomputed", learning_rate='auto')
        tsne_fit = tsne.fit(dist)
        dt = pd.DataFrame(tsne_fit.embedding_, columns=['PC1', 'PC2', 'PC3'])
        dt['sample'] = sample_names
        axis_names = ["t-SNE 1", "t-SNE 2", "t-SNE 3"]

    # generate plot dataframe
    df = dt.join(merged_metadata.set_index('sample'), on='sample', how='inner')
    df.index = df.index.map(str) # this step is important
    df.fillna('NA', inplace=True)

    # reorder df
    ess_cols = ['sample','collection_date','collection_site']
    df = df[ess_cols + [i for i in df.columns if i not in ess_cols]]
    # write out ordination_result
    df.to_csv("{0}/{1}_ordination_result.txt".format(outpath, analysis), index=False, sep="\t")

    logger.info('Plotting')
    hovertmp = generate_hovertemplate(list(df.columns), method=analysis)
    if column == "default":
        # plot 3d
        default_plot(df, analysis, intersect_len, outpath, axis_names, hovertmp, voc=voc, mode="3d")
        # plot 2d
        default_plot(df, analysis, intersect_len, outpath, axis_names, hovertmp, voc=voc, mode="2d")
    else:
        # others as column
        # plot 3d
        other_col_plot(df, column, analysis, intersect_len, outpath, axis_names, hovertmp, voc=voc, mode="3d")
        # plot 2d
        other_col_plot(df, column, analysis, intersect_len, outpath, axis_names, hovertmp, voc=voc, mode="2d")
    
    logger.info('Plotting is done!')

if __name__ == '__main__':
    cli()
