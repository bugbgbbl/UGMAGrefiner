# under python 3 

import argparse
#import subprocess
#from calendar import c
#from hashlib import new
import logging
#from platform import node
import time
from Bio import SeqIO
import re
from multiprocessing import Pool
import psutil
import os

#from numpy import cov, real_if_close
#from bidirectionalmap.bidirectionalmap import BidirectionalMap
import sys
#from collections import defaultdict
from igraph import *
import csv
import itertools as it
from typing import List, Dict
from sklearn.mixture import GaussianMixture
#from sklearn.cluster import DBSCAN
import numpy as np

ap = argparse.ArgumentParser(description="""UGMAGrefiner is a tool which 1) recruit unitigs to MAGs to improve compeleteness of MAGs which assembled by metaSPAdes and binning by existing binning tools, it is able to assign edges to multiple bins. 2) identify whether there are multiple genomes mixed in one MAG and get the unitigs of the different part of these genomes. 
UGMAGrefiner uses the connectivity and coverage information from assembly graphs to improve completeness of MAGs from exsiting binning tools on unitig level and to identify different part of similar genomes mixed in one MAG.""")

ap.add_argument("--edges", required=True, help="path to the assembly_graph.fastg file")
ap.add_argument("--graph", required=True, help="path to the assembly_graph_after_simplification.gfa file")
ap.add_argument("--paths", required=True, help="path to the contigs.paths file")
ap.add_argument("--binned", required=True, help="path to the .csv file with the initial binning output from an existing tool")
ap.add_argument("--output", required=True, help="path to the output folder")
ap.add_argument("--prefix", required=False, default='', help="prefix for the output file")
ap.add_argument("--depth", required=False, type=int, default=5, help="maximum depth for the breadth-first-search. [default: 5]")
ap.add_argument("--threshold", required=False, type=float, default=1.5, help="threshold for determining inconsistent vertices. [default: 1.5]")
ap.add_argument("--delimiter", required=False, type=str, default=",", help="delimiter for input/output results [default: , (comma)]")
ap.add_argument("--nthreads", required=False, type=int, default=8, help="number of threads to use. [default: 8]")
ap.add_argument("--gfaout", action="store_true", help="whether output gfa file for each MAG")
#ap.add_argument("--trueContigSource", required=False, type=str, help="ture contig source used to calculate the accuracy of result")

args = vars(ap.parse_args())

edges_file = args["edges"]
#edges_file="/home/xiangbaoyu/project/pwsT2D/7-assembly/metaspades/time2/assembly_graph.fastg"
assembly_graph_file = args["graph"]
#assembly_graph_file="/home/xiangbaoyu/project/pwsT2D/7-assembly/metaspades/time2/assembly_graph_after_simplification.gfa"
contig_paths = args["paths"]
#contig_paths="/home/xiangbaoyu/project/pwsT2D/7-assembly/metaspades/time2/contigs.paths"
contig_bins_file = args["binned"]
#contig_bins_file="/home/xiangbaoyu/project/pwsT2D/21-minerefine/bininput/time2_refine.csv"
output_path = args["output"]
#output_path="./"
prefix = args["prefix"]
#prefix="test"
depth = args["depth"]
#depth=6
threshold = args["threshold"]
#threshold=1.5
delimiter = args["delimiter"]
#delimiter=","
nthreads = args["nthreads"]
#nthreads=8
#true_contigSource = args["trueContigSource"]
gfaout = args["gfaout"]

# edges_file="/home/xiangbaoyu/project/pwsT2D/28-GD02/1-metaspades/ERR1753684/assembly_graph.fastg"
# assembly_graph_file="/home/xiangbaoyu/project/pwsT2D/28-GD02/1-metaspades/ERR1753684/assembly_graph_after_simplification.gfa"
# contig_paths="/home/xiangbaoyu/project/pwsT2D/28-GD02/1-metaspades/ERR1753684/contigs.paths"
# contig_bins_file="/home/xiangbaoyu/project/pwsT2D/28-GD02/3-minerefine/test/ERR1753684_initial_contig_bins.csv"
# output_path="./"
# prefix="test"
# depth=6
# threshold=1.5
# delimiter=","
# nthreads=16

# Setup timer and mem
start_time = time.time()
pid = os.getpid()
p = psutil.Process(pid)
info_start = p.memory_full_info().uss/1048576

# Setup logger
#-----------------------
logger = logging.getLogger('UGMAGrefiner')
logger.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
consoleHeader = logging.StreamHandler()
consoleHeader.setFormatter(formatter)
consoleHeader.setLevel(logging.INFO)
logger.addHandler(consoleHeader)

# Setup output path for log file
#---------------------------------------------------

fileHandler = logging.FileHandler(prefix+"UGMAGrefiner.log")
fileHandler.setLevel(logging.DEBUG)
fileHandler.setFormatter(formatter)
logger.addHandler(fileHandler)


logger.info("Welcome to UGMAGrefiner: Refined and identify genome specific region from Binning result of Metagenomic sequencing using unitig level Assembly Graphs.")
logger.info("This version of UGMAGrefiner makes use of the assembly graph produced by metaSPAdes which is based on the de Bruijn graph approach.")

logger.info("Input arguments:")
logger.info("edges file: "+edges_file)
logger.info("Assembly graph file: "+assembly_graph_file)
logger.info("Contig paths file: "+contig_paths)
logger.info("Existing binning output file: "+contig_bins_file)
logger.info("Final binning output file: "+output_path)
logger.info("Prefix: "+ prefix)
logger.info("Depth: "+str(depth))
logger.info("Threshold: "+str(threshold))
logger.info("Number of threads: "+str(nthreads))

logger.info("UGMAGrefiner started")

start_time = time.time()


def get_edge_info(edges_file:str):
    '''
    # Get length, coverage and translate dic of edges order to edges num

    edge_lengths: dic[str, int]
    edge_coverages: dic[str, float]
    edge_name: dic[str, str]
    '''

    edge_lengths = {}
    edge_coverages = {}
    edge_name={}

    try:
        for index, record in enumerate(SeqIO.parse(edges_file, "fasta")):
            start = 'EDGE_'
            end = '_length'
            edge_num = re.search('%s(.*?)%s' % (start, end), record.id).group(1)
            if(edge_num not in edge_name):
                edge_name[edge_num] = str(re.search('^(.*?)[\'\:\;]',record.id).group(1))

                start = '_length_'
                end = '_cov'
                length = int(re.search('%s(.*?)%s' % (start, end), record.id).group(1))
                
                start = '_cov_'
                end = "[';:]"
                coverage = float(re.search('%s(.*?)%s' % (start, end), record.id).group(1))
                
                edge_lengths[edge_num] = length
                edge_coverages[edge_num] = coverage
        logger.info("Total number of edges: "+str(len(edge_lengths)))
        return [edge_lengths, edge_coverages,edge_name]

    except:
        logger.error("Please make sure that the correct path to the assembly_graph.fastg file is provided.")
        logger.info("Exiting UGMAGrefiner... Bye...!")
        sys.exit(1)

def get_cotigs_info(contig_paths:str):
    '''
    # Get contig length, coverage, and paths from contigs.paths  

    paths: dic[str, list[str]]
    contig_lengths: dic[str, int]
    contig_coverages: dic[str, float]
    edge_contigs: dic[str, str]
    '''

    contig_lengths = {}
    contig_coverages = {}
    #record contigs' composition
    paths = {}
    #record edge's source
    edge_contigs = {}
    #contigs' name
    #contig_names = {}

    current_contig_num = ""

    try:
        with open(contig_paths) as file:
            name = file.readline()
            path = file.readline()
            
            while name != "" and path != "":
                #deal with ; in path line    
                while ";" in path:
                    path = path[:-2]+","+file.readline()
                
                #get contigs' number, length, coverage from it's name, only can be used in spades output
                start = 'NODE_'
                end = '_length_'
                contig_num = re.search('%s(.*)%s' % (start, end), name).group(1)

                edges = path.rstrip().split(",")

                if current_contig_num != contig_num:
                    start = '_length_'
                    end = '_cov'
                    length = int(re.search('%s(.*)%s' % (start, end), name).group(1))
        
                    start = '_cov_'
                    end = ''
                    coverage = float(re.search('%s(.*)%s' % (start, end), name).group(1))
                
                    contig_lengths[contig_num] = length
                    contig_coverages[contig_num] = coverage

                    #contig_names[contig_num] = name.strip()
                    current_contig_num = contig_num
                #record the contig's path
                if contig_num not in paths:
                    paths[contig_num] = [segment[:-1] for segment in edges]
                #record each edge's contigs source
                for edge in edges:
                    if edge not in edge_contigs:
                        edge_contigs[edge] = set([contig_num])
                    else:
                        edge_contigs[edge].add(contig_num)
                name = file.readline()
                path = file.readline()
        logger.info("Total number of contigs available: "+str(len(contig_lengths)))       
        return [contig_lengths, contig_coverages, paths, edge_contigs]
        
    except:
        logger.error("Please make sure that the correct path to the contig paths file: contigs.paths is provided.")
        logger.info("Exiting UGMAGrefiner... Bye...!")
        sys.exit(1)

        
    

def construct_assembly_graph(assembly_graph_file:str,edge_name:Dict[str,str]):
    '''
    ## Construct the assembly graph

    assembly_graph: graph from igraph
    edge_name: dic[str, str]

    edge_node: dic[str, int]
    node_edge: dic[int, str]
    '''

    links = []
    edge_node = {}
    node_edge = {}
    #try:
        # Get links from assembly_graph_with_scaffolds.gfa
       
    # Create graph
    assembly_graph = Graph()

    numberOfedge=len(edge_name)
    # Add vertices
    assembly_graph.add_vertices(numberOfedge*2)
    # Name vertices
    j=0
    try:
        for i in edge_name:
            edge_node[i] = j
            node_edge[j] = i
            assembly_graph.vs[2*j]["id"]= str(j)
            assembly_graph.vs[2*j]["label"]= "start"
            assembly_graph.vs[2*j+1]["id"]= str(j)
            assembly_graph.vs[2*j+1]["label"]= "end"
            j = j + 1

        with open(assembly_graph_file) as file:
            line = file.readline()
            while line != "":
                # Identify lines with link information
                if "L" in line:
                    strings = line.split("\t")
                    #in link line of a gfa file, the first edge's end is connect to second edge's start, while "-" indicate its reverse complementary
                    # here, we use 2i indicate the start of edge, 2i + 1 indicate the end.  So we can record the start-to-start, start-to-end, end-to-end or end-to-start connection between two edges
                    if(strings[2]== "+"):
                        start = 2 * int(edge_node[strings[1]]) + 1
                    else:
                        start = 2 * int(edge_node[strings[1]])
                    if(strings[4]== "+"):
                        end = 2 * int(edge_node[strings[3]])
                    else:
                        end = 2 * int(edge_node[strings[3]]) + 1
                    links.append((start,end))
                line = file.readline()
        

        # Add edges to the graph
        assembly_graph.add_edges(links)
        assembly_graph.simplify(multiple=True, loops=False, combine_edges=None)
        logger.info("Total number of edges links in the assembly graph: "+str(len(links)))
        return assembly_graph,edge_node, node_edge

    except:
       logger.error("Please make sure that the correct path to the assembly graph file: assembly_graph_after_simplification.gfa is provided.")
       logger.info("Exiting UGMAGrefiner... Bye...!")
       sys.exit(1)
    
    
def get_bininput(contig_bins_file: str, paths: Dict[str, List[str]], contig_coverages: Dict[str, float],contig_lengths: Dict[str, int]):
    '''
    # Get  initial binning result

    n_bins: int
    bins_contigs: List[List[str]]
    bins_edges: List[List[str]]
    edges_bin_coverage: Dict[str, Dict[int, float]]  restore the mapping relationship between edges to bins and corresponding coverage.
    bins_average_coverage: list[int]
    '''
    try:
        all_bins_list = []

        with open(contig_bins_file) as csvfile:
            readCSV = csv.reader(csvfile, delimiter=delimiter)
            for row in readCSV:
                all_bins_list.append(row[1])
                
        bins_list = list(set(all_bins_list))
        bins_list.sort(key=int)

        n_bins = int(bins_list[-1]) #number of bins
        logger.info("Number of bins available in binning result: "+str(n_bins))

        bins_contigs = [[] for x in range(n_bins)]
        bins_edges = [[] for x in range(n_bins)]
        edges_bin_coverage = {}  #Dict[str, Dict[int, float]]
        with open(contig_bins_file) as contig_bins:
            readCSV = csv.reader(contig_bins, delimiter=delimiter)
            for row in readCSV:
                start = 'NODE_'
                end = '_length_'
                contig_num = re.search('%s(.*?)%s' % (start, end), row[0]).group(1)
                
                bin_num = int(row[1])-1
                bins_contigs[bin_num].append(contig_num)
                bins_edges[bin_num].extend(paths[contig_num])
                #for edges in each contigs, if it has aleardy been added to edge_bin_coverage by other contigs, sum it's coverage, if not, record it.
                for cedge in paths[contig_num]:   
                    if cedge not in edges_bin_coverage:
                        edges_bin_coverage[cedge] = {}
                    if(bin_num not in edges_bin_coverage[cedge]):
                        edges_bin_coverage[cedge][bin_num] = contig_coverages[contig_num]
                    else:
                        edges_bin_coverage[cedge][bin_num] += contig_coverages[contig_num] 
            for i in range(len(bins_edges)):
                bins_edges[i] = list(set(bins_edges[i]))

        bins_contigs_len_sum=[0 for x in range(n_bins)]
        bins_contigs_lengthcoverage=[0 for x in range(n_bins)] 
        bins_average_coverage = [0 for x in range(n_bins)]  
        for bin in range(len(bins_contigs)):
            for contig in bins_contigs[bin]:
                bins_contigs_len_sum[bin] += contig_lengths[contig]
                bins_contigs_lengthcoverage[bin] += contig_lengths[contig] * contig_coverages[contig]
            bins_average_coverage[bin] = bins_contigs_lengthcoverage[bin] / bins_contigs_len_sum[bin]
        return [n_bins, bins_contigs, bins_edges,edges_bin_coverage,bins_average_coverage]
        
    except:
        logger.error("Please make sure that the correct path to the binning result file is provided and it is having the correct format")
        logger.info("Exiting UGMAGrefiner... Bye...!")
        sys.exit(1)    
    



def get_binned_and_unbinned_edges(n_bins: int, bins_edges: List[List[str]],edge_lengths:Dict[str,int]):
    '''
    # Get binned and unbinned edges

    n_bins: number of bins
    bins_edges: List[List[str]] edges(num) in each bins(bins name)
    edge_lengths: Dict[str, int]

    binned_edges: List[str]
    unbinned_edges: List[str]
    '''

    binned_edges = []

    for n in range(n_bins):
        binned_edges = binned_edges+[edge for edge in bins_edges[n]]
    binned_edges=list(set(binned_edges))
    unbinned_edges = list(set(edge_lengths).difference(set(binned_edges)))

    binned_edges.sort()
    unbinned_edges.sort()

    logger.info("Number of binned edges: "+str(len(binned_edges)))
    logger.info("Total number of unbinned edges: "+str(len(unbinned_edges)))

    return [binned_edges, unbinned_edges]




def runBFS(edge:str, node: int, edge_coverages: Dict[str,float], edges_bin_coverage: Dict[str,Dict[int,float]], binned_edges: List[str],node_edge: Dict[int,str],assembly_graph,bins_average_coverage: List[int], threhold=depth,):
    '''
    # The BFS function to search labelled nodes
    for a given node in graph, find its all neighbors who have been labelled in a bredth first way under a given depth

    return: set((edge:int, activate_edge:int, bin_num:int, depth:int, activate_edge_coverage/edge_coverage:float))    
    '''
    if(edge_coverages[edge]==0):
        logger.info("Something wrong with node coverage(equal to 0): "+str(node))
        return 0
    queue = []
    visited = set()
    queue.append(node)
    depth = {}
    
    depth[node] = 0
    
    labelled_nodes = set()
    
    while (len(queue) > 0):
        active_node = queue.pop(0)
        visited.add(active_node)

        active_node_num=int(active_node)//2
        active_edge = node_edge[active_node_num]
        #check if active node in binned_edges, if so, record it in labelled_nodes along with bins, depth, and foldchange of coverage
        if active_edge in binned_edges and len(visited) > 1:

            for i in edges_bin_coverage[active_edge]:
                labelled_nodes.add((int(node_edge[node//2]), int(node_edge[active_node//2]), i, 
                depth[active_node], bins_average_coverage[i] / edge_coverages[node_edge[node//2]]))

        #for node not binned, search its neighbor from its other side
        # elif((len(visited) == 1) or (len(assembly_graph.neighbors(active_node, mode="ALL")) == 1 and edge_coverages[active_edge] < 1.5*edge_coverages[edge]) or (edge_lengths[active_edge]< 150 and edge_coverages[active_edge] > edge_coverages[edge])):
        elif((len(visited) == 1) or (len(assembly_graph.neighbors(active_node, mode="ALL")) == 1 )):
            is_end=int(active_node)%2
            if(is_end==1):
                to_neighbor_num=int(active_node)-1
            else:
                to_neighbor_num=int(active_node)+1
            for neighbour in assembly_graph.neighbors(to_neighbor_num, mode="ALL"):
                if neighbour not in visited:
                    depth[neighbour] = depth[active_node] + 1
                    if depth[neighbour] > threhold:
                        continue
                    queue.append(neighbour)
    return labelled_nodes

def runBFS_for_newstrain(edge:str, node: int, edge_coverages: Dict[str,float], new_strain_edges_bin_cluster,node_edge: Dict[int,str],assembly_graph,threhold=depth,):
    '''
    # The BFS function to search new_strain nodes
    for a given node in graph, find its all neighbors who have been in new_strain in a bredth first way under a given depth

    return: set((edge:int, activate_edge:int, bin_num:int, depth:int, activate_edge_coverage/edge_coverage:float))    
    '''
    if(edge_coverages[edge]==0):
        logger.info("Something wrong with node coverage(equal to 0): "+str(node))
        return 0
    queue = []
    visited = set()
    queue.append(node)
    depth = {}
    
    depth[node] = 0
    labelled_nodes = set()

    while (len(queue) > 0):
        active_node = queue.pop(0)
        visited.add(active_node)

        active_node_num=int(active_node)//2
        active_edge = node_edge[active_node_num]
        #check if active node in new_strain_edges, if so, record it in labelled_nodes along with bins, depth, and foldchange of coverage
        if active_edge in new_strain_edges_bin_cluster and len(visited) > 1:
            labelled_nodes.add((int(node_edge[node//2]), int(node_edge[active_node//2]), str(new_strain_edges_bin_cluster[active_edge][0])+"_"+str(new_strain_edges_bin_cluster[active_edge][1]), depth[active_node], new_strain_edges_bin_cluster[active_edge][2] / edge_coverages[node_edge[node//2]]))

        #for node not binned, search its neighbor from its other side
        elif((len(visited) == 1) or (len(assembly_graph.neighbors(active_node, mode="ALL")) == 1)):
            is_end=int(active_node)%2
            if(is_end==1):
                to_neighbor_num=int(active_node)-1
            else:
                to_neighbor_num=int(active_node)+1
            for neighbour in assembly_graph.neighbors(to_neighbor_num, mode="ALL"):
                if neighbour not in visited:
                    depth[neighbour] = depth[active_node] + 1
                    if depth[neighbour] > threhold:
                        continue
                    queue.append(neighbour)
    return labelled_nodes

def merge_left_right(left_labelled_node_list,right_labelled_node_list,num_of_loop):
    '''
    # merge two labeled_nodes from one node' head and tail
    for two given labeled_nodes, merge them by calculating the sum coverage of each MAG and maintain the bigger one if the MAG exist in both side, if the MAG only shows in one side, mantain it in the result.

    return:
    source_bins: like labeled_nodes, is a set as: set((edge:int, activate_edge:int, bin_num:int, depth:int, activate_edge_coverage/edge_coverage:float))  
    source_bins_depth: a dict contain the depth for each MAG
    '''
    #store thesource bin's contribute for each edge  
    l_source_bins={}
    r_source_bins={}
    source_bins={}
    l_source_bins_depth={}
    r_source_bins_depth={}
    source_bins_depth={}

    #merge left and right labelled node
    if(left_labelled_node_list):
        tmp_edge=set()
        for i in left_labelled_node_list:
            if(i[1] not in tmp_edge):
                if(i[2] not in l_source_bins):
                    l_source_bins[i[2]]=0
                    l_source_bins_depth[i[2]]=10000
                tmp_edge.add(i[1])
                #calculate each connection's contribute to this edge. contribution = each edge's coverage
                l_source_bins[i[2]] += edge_coverages[str(i[1])]
                c_depth = (i[3]+(num_of_loop -1)*depth)
                if(c_depth < l_source_bins_depth[i[2]]):
                    l_source_bins_depth[i[2]] = c_depth
    if(right_labelled_node_list):
        tmp_edge=set()
        for i in right_labelled_node_list:
            if(i[1] not in tmp_edge):
                if(i[2] not in r_source_bins):
                    r_source_bins[i[2]]=0
                    r_source_bins_depth[i[2]]=10000
                tmp_edge.add(i[1])
                #calculate each connection's contribute to this edge. contribution = each edge's coverage
                r_source_bins[i[2]] += edge_coverages[str(i[1])]
                c_depth = (i[3]+(num_of_loop -1)*depth)
                if(c_depth < r_source_bins_depth[i[2]]):
                    r_source_bins_depth[i[2]] = c_depth
    #for bins in each side, obtain the max coverage as the neighbor bins' coverage result
    for i in l_source_bins.keys():
        if((i in r_source_bins and l_source_bins[i] > r_source_bins[i]) or (i not in r_source_bins) ):
            r_source_bins[i]=l_source_bins[i]
            r_source_bins_depth[i]=l_source_bins_depth[i]
    source_bins = r_source_bins
    source_bins_depth = r_source_bins_depth
    return source_bins,source_bins_depth

def givelabel(tolabel_edge: str):
    '''
    # get label from runBFS result
    
    input: tolabel_edge

    output: a list contain: [tolabel_edge, [current_bin], coverage_foldchange, [source_bins_depth[current_bin]],source_bins] 
    '''
    left_labelled_node_list=[item for item in runBFS(tolabel_edge, 2*edge_node[tolabel_edge], edge_coverages, edges_bin_coverage, binned_edges,node_edge,assembly_graph,bins_average_coverage,threhold=depth)]
    right_labelled_node_list=[item for item in runBFS(tolabel_edge, 2*edge_node[tolabel_edge] + 1, edge_coverages, edges_bin_coverage, binned_edges,node_edge,assembly_graph,bins_average_coverage,threhold=depth)]

    source_bins, source_bins_depth = merge_left_right(left_labelled_node_list, right_labelled_node_list,num_of_loop)

    from_label_number = len(source_bins)   
    if(from_label_number == 1):
        current_bin = list(source_bins.keys())[0]
        coverage_foldchange = edge_coverages[tolabel_edge] / bins_average_coverage[current_bin]
        if(edge_coverages[tolabel_edge]>1.5*source_bins[current_bin]):
            return [tolabel_edge, [-1], coverage_foldchange,[source_bins_depth[current_bin]],source_bins]
        #elif(0.1*bins_average_coverage[current_bin]<=edge_coverages[tolabel_edge]<=1.5*source_bins[current_bin]):
        elif(edge_coverages[tolabel_edge]<=1.5*source_bins[current_bin]):
            return [tolabel_edge, [current_bin], coverage_foldchange, [source_bins_depth[current_bin]],source_bins]
        
    elif(from_label_number > 1):
        bin_combinations = []
        for i in range(len(source_bins)):
            bin_combinations += list(it.combinations(source_bins, i+1))
        min_diff = sys.maxsize
        min_diff_combination = -1 
        min_comb_cov_total = 0
        min_cov_of_edges = 0
        for combination in bin_combinations:

            comb_cov_total = 0
            cov_of_edges = 0

            for i in range(len(combination)):
                comb_cov_total += bins_average_coverage[combination[i]]
                cov_of_edges += source_bins[combination[i]]
            cov_diff = min(abs(cov_of_edges-edge_coverages[tolabel_edge]),abs(comb_cov_total-edge_coverages[tolabel_edge]))

            if cov_diff < min_diff:
                min_diff = cov_diff
                min_diff_combination = combination
                min_comb_cov_total = comb_cov_total
                min_cov_of_edges = cov_of_edges

        if(min_diff_combination!=-1):
            foldchange_tmp=edge_coverages[tolabel_edge]/min_comb_cov_total
            if(len(min_diff_combination) == 1 and edge_coverages[tolabel_edge] < 1.5*min_cov_of_edges):
                return [tolabel_edge, min_diff_combination, foldchange_tmp, [source_bins_depth[x] for x in min_diff_combination],source_bins]
            elif(len(min_diff_combination) > 1 and foldchange_tmp < 1.5):           
                return [tolabel_edge, min_diff_combination, foldchange_tmp, [source_bins_depth[x] for x in min_diff_combination],source_bins]
            else:
                return [tolabel_edge, [-1], foldchange_tmp, [source_bins_depth[x] for x in min_diff_combination],source_bins]
                            

def extend_cluster(tolabel_edge: str):
    '''
    # extend identified clusters to calculate extend rate

    input: tolabel_edge from identified clusters

    output: a list like: [tolabel_edge, min_diff_combination, foldchange_tmp, [source_bins_depth[x] for x in min_diff_combination]]
    '''
    left_labelled_node_list=[item for item in runBFS_for_newstrain(tolabel_edge, 2*edge_node[tolabel_edge], edge_coverages, new_strain_edges_bin_cluster, node_edge,assembly_graph,threhold=5)]
    right_labelled_node_list=[item for item in runBFS_for_newstrain(tolabel_edge, 2*edge_node[tolabel_edge] + 1, edge_coverages, new_strain_edges_bin_cluster, node_edge,assembly_graph,threhold=5)]
    #return set: [to_edge, source_edge, bin"_"clsuter, depth, cluster_mean_coverage/edge_coverage]
    from_label_number= 0 

    #store thesource bin's contribute for each edge  
    l_source_bins={}
    r_source_bins={}
    source_bins={}
    l_source_bins_depth={}
    r_source_bins_depth={}
    source_bins_depth={}

    #merge left and right labelled node
    if(left_labelled_node_list):
        tmp_edge=set()
        for i in left_labelled_node_list:
            if(i[1] not in tmp_edge):
                if(i[2] not in l_source_bins):
                    l_source_bins[i[2]]=0
                    l_source_bins_depth[i[2]]=10000
                tmp_edge.add(i[1])
                #calculate each connection's contribute to this edge. contribution = each edge's coverage
                l_source_bins[i[2]] += edge_coverages[str(i[1])]
                c_depth = (i[3]+(num_of_loop -1)*depth)
                if(c_depth < l_source_bins_depth[i[2]]):
                    l_source_bins_depth[i[2]] = c_depth
    if(right_labelled_node_list):
        tmp_edge=set()
        for i in right_labelled_node_list:
            if(i[1] not in tmp_edge):
                if(i[2] not in r_source_bins):
                    r_source_bins[i[2]]=0
                    r_source_bins_depth[i[2]]=10000
                tmp_edge.add(i[1])
                #calculate each connection's contribute to this edge. contribution = each edge's coverage
                r_source_bins[i[2]] += edge_coverages[str(i[1])]
                c_depth = (i[3]+(num_of_loop -1)*depth)
                if(c_depth < r_source_bins_depth[i[2]]):
                    r_source_bins_depth[i[2]] = c_depth
    #only obtain the bin that occur in both side 
    for i in list(l_source_bins.keys()):
        if((i in r_source_bins and l_source_bins[i] > r_source_bins[i]) or (i not in r_source_bins) ):
            r_source_bins[i]=l_source_bins[i]
            r_source_bins_depth[i]=l_source_bins_depth[i]
    source_bins = r_source_bins
    source_bins_depth = r_source_bins_depth


    from_label_number = len(source_bins)   
    if(from_label_number == 1):
        current_bin_cluster=list(source_bins.keys())[0]
        current_bin = int(list(source_bins.keys())[0].split("_")[0])
        current_cluster = int(list(source_bins.keys())[0].split("_")[1])
        coverage_foldchange = edge_coverages[tolabel_edge] / new_strain_foldchange[current_bin][current_cluster]
        if(edge_coverages[tolabel_edge]>2*source_bins[current_bin_cluster]):
            return [tolabel_edge, [-1], coverage_foldchange,[source_bins_depth[current_bin_cluster]]]
        elif(0.1*new_strain_foldchange[current_bin][current_cluster]<=edge_coverages[tolabel_edge]<=2*source_bins[current_bin_cluster]):
            return [tolabel_edge, [list(source_bins.keys())[0]], coverage_foldchange, [source_bins_depth[current_bin_cluster]]]
        
    elif(from_label_number > 1):
        bin_combinations = []
        for i in range(len(source_bins)):
            bin_combinations += list(it.combinations(source_bins, i+1))
        min_diff = sys.maxsize
        min_diff_combination = -1 
        min_comb_cov_total = 0
        min_cov_of_edges = 0
        for combination in bin_combinations:

            comb_cov_total = 0
            cov_of_edges = 0

            for i in range(len(combination)):
                current_bin = int(combination[i].split("_")[0])
                current_cluster = int(combination[i].split("_")[1])
                comb_cov_total += new_strain_foldchange[current_bin][current_cluster]
                cov_of_edges += source_bins[combination[i]]
            cov_diff = min(abs(cov_of_edges-edge_coverages[tolabel_edge]),abs(comb_cov_total-edge_coverages[tolabel_edge]))

            if cov_diff < min_diff:
                min_diff = cov_diff
                min_diff_combination = combination
                min_comb_cov_total = comb_cov_total
                min_cov_of_edges = cov_of_edges

        if(min_diff_combination!=-1):
            foldchange_tmp=edge_coverages[tolabel_edge]/min_comb_cov_total
            if(len(min_diff_combination) == 1 and 0.1*min_comb_cov_total < edge_coverages[tolabel_edge] < 2*min_cov_of_edges):
                return [tolabel_edge, min_diff_combination, foldchange_tmp, [source_bins_depth[x] for x in min_diff_combination]]
            if(len(min_diff_combination) > 1 and 0.1 < foldchange_tmp < 2):           
                return [tolabel_edge, min_diff_combination, foldchange_tmp, [source_bins_depth[x] for x in min_diff_combination]]
            elif((len(min_diff_combination) == 1 and edge_coverages[tolabel_edge] >= 2*min_cov_of_edges) or (len(min_diff_combination) > 1 and foldchange_tmp >= 2)):
                #high_coverage_edges.add(tolabel_edge)
                return [tolabel_edge, [-1], foldchange_tmp, [source_bins_depth[x] for x in min_diff_combination]]


def write_result(output_path, new_add_edges, edge_name):
    output_bins=[]
    #remove new add edges > 1.5Mb, which may be a new genome instead of a part of the original bins
    for bin in new_add_edges:
        if(len_sum(new_add_edges[bin],edge_length_info= edge_lengths)< 1500000):
            for edge in new_add_edges[bin]:
                output_bins.append([edge_name[str(edge)], int(bin)+1])
    output_file = output_path + prefix + 'UGMAGrefiner_output.csv'
    with open(output_file, mode='w') as output_file:
        output_writer = csv.writer(output_file, delimiter=delimiter, quotechar='"', quoting=csv.QUOTE_MINIMAL)
    
        for row in output_bins:
            output_writer.writerow(row)

    logger.info("Final binning results can be found at "+str(output_file.name))


def len_sum(edge_list, edge_length_info):
    '''
    # calculate the total length of edges exist in one edge_list
    '''
    length_sum = 0
    for edge in edge_list:
        length_sum = length_sum + edge_length_info[str(edge)]
    return length_sum


def check_labeled_nodes(binned_edges: List[str], bins_average_coverage: List[float],edge_coverages,edges_bin_coverage:Dict[str, Dict[int, float]],unbinned_edges:List[str]):
    '''
    # check binned nodes, remove the label(s) of nodes with abnormal high coverage 

    return: 
    binned_edges: new binned_edges with some nodes removed
    removed_edge: List[[str,[int]]], [[edge,[bin1,bin2,.....]]....]
    '''
    removed_edge = []
    #find edges with edge depth more than 1.5 times of bins depth
    to_refine_edges=set()
    for edge in binned_edges:
        coverage_sum=0
        for bin in edges_bin_coverage[edge]:
            coverage_sum+=bins_average_coverage[bin]
        if(coverage_sum > 0 and (edge_coverages[edge] / coverage_sum) > 1.5 ):
            to_refine_edges.add(edge)
            unbinned_edges.append(edge)
            removed_edge.append([edge,list(edges_bin_coverage[edge].keys())])                
            #edges_bin_coverage.pop(edge)   
    for key in to_refine_edges:
        edges_bin_coverage.pop(key) 

    binned_edges=list(filter(lambda x : x not in to_refine_edges, binned_edges))
    return binned_edges, removed_edge


def dogmm(new_add_edges:Dict[int,List[str]]):
    '''
    # do GaussianMixture based on unitigs' coverage for each MAG

    return:
    new_strain_bin_cluster_edges, Dict[int, Dict[int, Dict[str, List[float, int, int]]]]
    new_strain_edges_bin_cluster, Dict[str, List[bin, cluster, gmmmeans]]
    new_strain_length, Dict[str, int]]
    new_strain_foldchange, Dict[int, Dict[str, int]]
    multiple_source_edges, Dict[edge:[bin_cluster]]
    '''
    gmm_bin_category_edge={}
    new_strain_bin_cluster_edges={} #Dict[int, Dict[int, Dict[str, List[float, int, int]]]]
    new_strain_edges_bin_cluster={} #Dict[str, List[bin, cluster, gmmmeans]]
    new_strain_length={} #Dict[str, int]]
    new_strain_foldchange={} #Dict[int, Dict[str, int]]
    multiple_source_edges = {} #Dict[edge:[bin_cluster]]

    for bin in new_add_edges:
        gmm_bin_category_edge[bin]={}
        new_strain_bin_cluster_edges[bin]={}
        new_strain_foldchange[bin]={}
        to_gmm_list_foldchange =[]
        to_gmm_list_name=[]
        #only do gmm for edges with length > 100
        for edge in new_add_edges[bin]:
            if(edge_lengths[edge] >100):
                to_gmm_list_foldchange.extend([edge_coverages[edge]])
                to_gmm_list_name.append(edge)
        to_gmm_array = np.array(to_gmm_list_foldchange,dtype=float)
        to_gmm_array = to_gmm_array.reshape(-1,1)
        #only do gmm when number of edge > 20
        if(len(to_gmm_array) > 20):
            n_components = np.arange(2,min(10,len(to_gmm_array)))
            models = [GaussianMixture(n, random_state=0).fit(to_gmm_array) for n in n_components]
            models_bic = [m.bic(to_gmm_array) for m in models]
            models_aic = [m.aic(to_gmm_array) for m in models]
            n_to_gmm = models_bic.index(min(models_bic)) + 2

            gmm = GaussianMixture(n_components=n_to_gmm)
            gmm_model = gmm.fit(to_gmm_array) 
            gmm_predict = gmm.predict(to_gmm_array)
            # for i in range(len(gmm.means_)):
            #     gmm_out_file.write("gmmmean:%d,%d,%f\n" %(bin+1,i,gmm.means_[i][0]))

            merge_cluster = {} #Dict[newcluster:oldcluster]
            for i in range(len(gmm_predict)):
                if(gmm_predict[i] not in gmm_bin_category_edge[bin]):
                    gmm_bin_category_edge[bin][gmm_predict[i]]=[]
                gmm_bin_category_edge[bin][gmm_predict[i]].append(to_gmm_list_name[i])
                #if new cluster's foldchange are in the range (0.7,1.5) of old cluster, merge these two clusters
                if(gmm_predict[i] in merge_cluster):
                    c_cluster = merge_cluster[gmm_predict[i]]
                else:
                    c_cluster = gmm_predict[i]
                    if(c_cluster not in new_strain_foldchange[bin]):
                        for cluster in new_strain_foldchange[bin]:
                            if( 0.7 < new_strain_foldchange[bin][cluster] / gmm.means_[c_cluster][0] < 1.5):
                                merge_cluster[c_cluster] = cluster
                                c_cluster = cluster
                                break
                        if(c_cluster not in new_strain_bin_cluster_edges[bin]):
                            new_strain_bin_cluster_edges[bin][c_cluster] = {}
                            new_strain_length[str(bin)+"_"+str(c_cluster)] = 0
                            new_strain_foldchange[bin][c_cluster] = gmm.means_[gmm_predict[i]][0]

                new_strain_bin_cluster_edges[bin][c_cluster][to_gmm_list_name[i]] =[to_gmm_list_foldchange[i], edge_lengths[to_gmm_list_name[i]], gmm.means_[c_cluster][0]]
                if(to_gmm_list_name[i] not in new_strain_edges_bin_cluster):
                    new_strain_edges_bin_cluster[to_gmm_list_name[i]]=[bin, c_cluster, gmm.means_[c_cluster][0]]
                else:
                    if(to_gmm_list_name[i] not in multiple_source_edges):
                        multiple_source_edges[to_gmm_list_name[i]] = [new_strain_edges_bin_cluster[to_gmm_list_name[i]]]
                        multiple_source_edges[to_gmm_list_name[i]].append([bin,c_cluster,gmm.means_[c_cluster][0]])
                    else:
                        multiple_source_edges[to_gmm_list_name[i]].append([bin,c_cluster,gmm.means_[c_cluster][0]])
                new_strain_length[str(bin)+"_"+str(c_cluster)] += edge_lengths[to_gmm_list_name[i]]
    #if there are edges in multiple cluster, put the edge in the longest cluster
    for edge in multiple_source_edges:
        tmp_length=0
        for bin_cluster in multiple_source_edges[edge]:
            if(tmp_length==0):
                max_genome=bin_cluster
                tmp_length=new_strain_length[str(bin_cluster[0])+"_"+str(bin_cluster[1])]
            elif(new_strain_length[str(bin_cluster[0])+"_"+str(bin_cluster[1])] >= tmp_length):
                tmp_length=new_strain_length[str(bin_cluster[0])+"_"+str(bin_cluster[1])]
                new_strain_length[str(max_genome[0])+"_"+str(max_genome[1])]=new_strain_length[str(max_genome[0])+"_"+str(max_genome[1])] - edge_lengths[edge]
                del new_strain_bin_cluster_edges[max_genome[0]][max_genome[1]][edge]
                max_genome=bin_cluster
            elif(new_strain_length[str(bin_cluster[0])+"_"+str(bin_cluster[1])] < tmp_length):
                new_strain_length[str(bin_cluster[0])+"_"+str(bin_cluster[1])]=new_strain_length[str(bin_cluster[0])+"_"+str(bin_cluster[1])] - edge_lengths[edge]
                del new_strain_bin_cluster_edges[bin_cluster[0]][bin_cluster[1]][edge]
        new_strain_edges_bin_cluster[edge]=[bin_cluster[0], bin_cluster[1], bin_cluster[2]]


    return new_strain_bin_cluster_edges, new_strain_edges_bin_cluster, new_strain_length, new_strain_foldchange,multiple_source_edges      

def getseq(seq_dict,outfile):
    '''
    # extract all newly identified unitigs in one fasta format file
    input: seq_dict: a dict contains the unitigs' name and sequences 
    '''
    with open(outfile,"w") as fout:
        with open(output_path + prefix + "newly_identified_cluster.csv","r") as new_strain_out_file:
            for line in new_strain_out_file:
                word = line.strip().split(",")
                if(word[2] in seq_dict):
                    fout.write(">"+word[2]+"\n"+seq_dict[word[2]]+"\n")

def getsubgraph(new_add_edges:Dict[int,List[str]],assembly_graph_file:str,bins_edges:List[List[str]]):
    '''
    # get newly added unitigs for each MAG in fasta format, and if gfaout is true, output gfa format file for all unitigs in each MAG 
    '''
    edge_seq = {}
    edge_link1={}
    edge_link2={}
    with open(assembly_graph_file,"r") as gfafile:
        for line in gfafile:
            word=line.strip().split("\t")
            if(line.startswith("S")):
                edge_seq[word[1]] = line
            elif(line.startswith("L")):
                edge_link1[word[1]] = [word[3],line]
                edge_link2[word[3]] = [word[1],line]
    for bin in new_add_edges:
        all_edges = set()
        all_edges.update(set(new_add_edges[bin]))
        with open(output_path + prefix + "_bin"+str(bin+1)+"added.fasta","w") as new_add_fastaout:
            for edge in all_edges:
                if(edge in edge_seq):
                    tmp_seq =  edge_seq[edge].strip().split("\t")[2]
                    new_add_fastaout.write(">%s\n%s\n" % (edge, tmp_seq))

        if(gfaout):
            with open(output_path + prefix + "_bin"+str(bin+1)+"added.gfa","w") as fout:
                all_edges.update(set(bins_edges[bin]))
                for edge in all_edges:
                    if(edge in edge_seq):
                        fout.write(edge_seq[edge])
                    if(edge in edge_link1 and edge_link1[edge][0] in all_edges):
                        fout.write(edge_link1[edge][1])
                    if(edge in edge_link2 and edge_link2[edge][0] in all_edges):
                        fout.write(edge_link2[edge][1])
def count_time(func):
    def int_time():
        start_time = time.time()
        func()
        over_time = time.time()
        total_time = over_time - start_time
        logger.info("running time: %s Second" % total_time)
    return int_time

def count_mem_info(func):
    def float_info():
        pid = os.getpid()
        p = psutil.Process(pid)
        info_start = p.memory_full_info().uss/1024
        func()
        info_end=p.memory_full_info().uss/1024
        logger.info("Memory used: "+str(info_end-info_start)+"KB")
    return float_info

##########initialize #######
if(not os.path.exists(output_path)):
    logger.info("output path not exist, creating")
    try:
        os.makedirs(output_path) 
    except:
        logger.info("cannot create output path, please check output_path")

new_add_edges = {} #Dict[int, List[str]],{bin:[edges]}
high_coverage_edges = set()
to_gmm_edge_bins={}
to_gmm_bin_edges = {}

edge_lengths, edge_coverages, edge_name = get_edge_info(edges_file)
contig_lengths, contig_coverages, paths, edge_contigs = get_cotigs_info(contig_paths)
assembly_graph, edge_node, node_edge = construct_assembly_graph(assembly_graph_file,edge_name)
n_bins, bins_contigs, bins_edges, edges_bin_coverage,bins_average_coverage = get_bininput(contig_bins_file, paths, contig_coverages,contig_lengths)
binned_edges, unbinned_edges = get_binned_and_unbinned_edges(n_bins, bins_edges,edge_lengths)
#remove label(s) of nodes with abnormal high coverages
binned_edges, refined_edge = check_labeled_nodes(binned_edges, bins_average_coverage,edge_coverages,edges_bin_coverage,unbinned_edges)


#loop to give label(s) to unlabeled nodes
#-----------------------------
previous_number_new_add_edges=0
number_new_add_edges=0
num_of_loop=0
to_gmm_edge_bins={}
while((number_new_add_edges == 0 or number_new_add_edges != previous_number_new_add_edges) and num_of_loop <= 10):

    c_new_add_edges=set()
    num_of_loop+=1

    #search each unbinned edges' label
    with Pool(nthreads) as my_pool:
        tmp_result = list(my_pool.imap_unordered(givelabel,unbinned_edges))
    
    #add result to new add edges
    for label_result in tmp_result:
        if(label_result and label_result[1][0] != -1):
            if(label_result[2] >= 0.7):
                if(len(label_result[1])>1):
                    #new_add_edge_info[label_result[0]]=""
                    for current_bin_index in range(len(label_result[1])):
                        current_bin = label_result[1][current_bin_index]
                        if(current_bin not in new_add_edges):
                            new_add_edges[current_bin]=[]
                        new_add_edges[current_bin].append(label_result[0])

                else:
                    if(label_result[1][0] not in new_add_edges):
                        new_add_edges[label_result[1][0]]=[]
                    new_add_edges[label_result[1][0]].append(label_result[0])
                    
            elif(num_of_loop == 1): #put all < 0.5 foldchange edge to a list, then used to do gmm
                if(len(label_result[4])>1):
                    for current_bin in list(label_result[4].keys()):
                        if(label_result[0] not in to_gmm_edge_bins):
                            to_gmm_edge_bins[label_result[0]]=set()
                        to_gmm_edge_bins[label_result[0]].add(current_bin)
                else:
                    if(label_result[0] not in to_gmm_edge_bins):
                        to_gmm_edge_bins[label_result[0]]=set()
                    to_gmm_edge_bins[label_result[0]].add(list(label_result[4].keys())[0])
    #update binned, unbinned, and edge_bin_coverages for next loop
    for cbin in new_add_edges:
        if(new_add_edges[cbin]):
            c_new_add_edges.update(set(new_add_edges[cbin]))
            binned_edges=list(set(binned_edges).union(set(new_add_edges[cbin])))
            unbinned_edges=list(set(unbinned_edges).difference(set(new_add_edges[cbin])))
            for cedge in new_add_edges[cbin]:
                if(cedge not in edges_bin_coverage):
                    edges_bin_coverage[cedge]={}
                edges_bin_coverage[cedge][cbin]=bins_average_coverage[cbin]

    #update this loop's number of new add edges
    previous_number_new_add_edges=number_new_add_edges
    number_new_add_edges=len(c_new_add_edges)
    
    logger.info("loop times:"+str(num_of_loop)+",previous new add edges"+str(previous_number_new_add_edges)+",now new add edgs:"+str(number_new_add_edges))
    logger.info("give label finished")
###############################################################
#remove edges in to_gmm_bin_edges which also in new_add_edges, save the result in a set to_gmm_edges
all_add_edges=set()
for cbin in new_add_edges:
    all_add_edges = all_add_edges.union(set(new_add_edges[cbin])) 
to_gmm_edges=set()
for edge in list(to_gmm_edge_bins.keys()):
    if(edge in all_add_edges):
        del to_gmm_edge_bins[edge]
    else:
        to_gmm_edges.add(edge)

#################################################
#get to_gmm_bin_edges from to_gmm_edge_bins
for edge in to_gmm_edge_bins:
    for bin in to_gmm_edge_bins[edge]:
        if(bin not in to_gmm_bin_edges):
            to_gmm_bin_edges[bin]=[]
        to_gmm_bin_edges[bin].append(edge)

###############################################
#cluster edges using GaussianMixture based on their coverage in each MAG
logger.info("clustering edges using GaussianMixture based on their coverage in each MAG") 
new_strain_bin_cluster_edges, new_strain_edges_bin_cluster, new_strain_length, new_strain_foldchange,multiple_source_edges = dogmm(to_gmm_bin_edges)
logger.info("clustering finished") 

######################################################
#extend each cluster to claculate extend rate
logger.info("extending each cluster to claculate extend rate") 
unbinned_edges=list(set(unbinned_edges).difference(set(new_strain_edges_bin_cluster.keys())))
with Pool(nthreads) as my_pool:
    tmp_result = list(my_pool.imap_unordered(extend_cluster,list(unbinned_edges)))

extend_bin_cluster_length={}
extend_bin_cluster_edges={}
for label_result in tmp_result:
    if(label_result and label_result[1][0] != -1):
        if(label_result[2] >= 0.5):
            if(len(label_result[1])>1):
                for current_bin_cluster_index in range(len(label_result[1])):
                    current_bin_cluster = label_result[1][current_bin_cluster_index]
                    if(current_bin_cluster not in extend_bin_cluster_length):
                        extend_bin_cluster_length[current_bin_cluster]=0
                        extend_bin_cluster_edges[current_bin_cluster]=set()
                    extend_bin_cluster_length[current_bin_cluster] += edge_lengths[label_result[0]]
                    extend_bin_cluster_edges[current_bin_cluster].add(label_result[0])
            else:
                if(label_result[1][0] not in extend_bin_cluster_length):
                    extend_bin_cluster_length[label_result[1][0]]=0
                    extend_bin_cluster_edges[label_result[1][0]]=set()
                extend_bin_cluster_length[label_result[1][0]] += edge_lengths[label_result[0]]
                extend_bin_cluster_edges[label_result[1][0]].add(label_result[0])

##############################
#write extend result
with open(output_path + prefix + "extend_bin_cluster.csv","w") as extend_file:
    for bin in new_strain_bin_cluster_edges:
        for group in new_strain_bin_cluster_edges[bin]:
            bin_cluster = str(bin)+'_'+str(group)            
            if(bin_cluster in extend_bin_cluster_length):
                if(new_strain_length[bin_cluster]+extend_bin_cluster_length[bin_cluster] > 10000):
                    extend_file.write("gmmextend:"+str(bin_cluster) + ","+ str(extend_bin_cluster_length[bin_cluster])+","+str(new_strain_length[bin_cluster])+"\n")
                    for edge in new_strain_bin_cluster_edges[bin][group]:
                        extend_file.write("%d,%d,%s,%f,%d,%f\n" % (bin+1,group, edge,new_strain_bin_cluster_edges[bin][group][edge][0],new_strain_bin_cluster_edges[bin][group][edge][1], new_strain_bin_cluster_edges[bin][group][edge][2]))
                    for edge in extend_bin_cluster_edges[bin_cluster]:
                        extend_file.write("%d,%d,%s,%f,%d,%f\n" % (bin+1,group, edge,edge_coverages[edge],edge_lengths[edge], new_strain_foldchange[bin][group]))
            if((bin_cluster not in extend_bin_cluster_length) and new_strain_length[bin_cluster] > 10000):
                extend_file.write("gmmextend:"+str(bin_cluster) + ","+ str(0)+","+str(new_strain_length[bin_cluster])+"\n")
                for edge in new_strain_bin_cluster_edges[bin][group]:
                    extend_file.write("%d,%d,%s,%f,%d,%f\n" % (bin+1,group, edge,new_strain_bin_cluster_edges[bin][group][edge][0],new_strain_bin_cluster_edges[bin][group][edge][1], new_strain_bin_cluster_edges[bin][group][edge][2]))

#get unitigs' sequences
seq_dict={} #Dict[str,str]
with open(assembly_graph_file,"r") as fgfa:
    for line in fgfa:
        if(line.startswith("S")):
            word = line.strip().split()
            seq_dict[word[1]]=word[2]        

#output newly identified clsuters
with open(output_path + prefix +"newly_identified_cluster.csv","w") as new_strain_out_file:
    for bin in new_strain_bin_cluster_edges:
        for group in new_strain_bin_cluster_edges[bin]:
            cbin_cluster = str(bin)+"_"+str(group)
            if(cbin_cluster in extend_bin_cluster_length and new_strain_length[cbin_cluster] > 0 ):
                extend_rate = extend_bin_cluster_length[cbin_cluster] / new_strain_length[cbin_cluster]
            else:
                extend_rate = 0
            cluster_type=""
            if(new_strain_length[cbin_cluster] > 1000000):
                cluster_type = "New" 
            elif(new_strain_length[cbin_cluster] > 100000 and new_strain_foldchange[bin][group] > 10):
                if(extend_rate > 1):
                    cluster_type="New"
                else:
                    cluster_type="Separated"
            if(cluster_type=="Separated" or cluster_type=="New"):
                for edge in new_strain_bin_cluster_edges[bin][group]:
                    new_strain_out_file.write("%d,%d,%s,%f,%d,%f,%s\n" % (bin+1,group, edge,new_strain_bin_cluster_edges[bin][group][edge][0],new_strain_bin_cluster_edges[bin][group][edge][1], new_strain_bin_cluster_edges[bin][group][edge][2], cluster_type))

                with open(output_path + prefix + "_"+str(bin+1) + "_" + str(group) + cluster_type+"_cluster.fasta","w") as newly_identified_cluster:
                    for edge in new_strain_bin_cluster_edges[bin][group]:
                        newly_identified_cluster.write(">%s\n%s\n" % (edge, seq_dict[edge]))

getseq(seq_dict, prefix+"_cluster.fasta")


getsubgraph(new_add_edges, assembly_graph_file, bins_edges)
write_result(output_path, new_add_edges, edge_name)

over_time = time.time()
total_time = over_time - start_time
logger.info("running time %d second" % total_time)
info_end=p.memory_full_info().uss/1048576
logger.info("used memory %d MB" % (info_end-info_start))
