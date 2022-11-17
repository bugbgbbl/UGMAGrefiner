
# UGMAGrefiner
Unitig level assembly Graph based Metagenome-assembled Genome refiner (UGMAGrefiner): a tool to increase completeness and resolution of metagenome-assembled genomes
# Dependencies
igraph
Bio
sklearn
numpy


# Install UGMAGrefiner
download UGMAGrefienr.py and all dependencies

# Usage
python UGMAGrefienr -h
usage: minerefine_SPAdes.py [-h] --edges EDGES --graph GRAPH --paths PATHS
                            --binned BINNED --output OUTPUT [--prefix PREFIX]
                            [--depth DEPTH] [--threshold THRESHOLD]
                            [--delimiter DELIMITER] [--nthreads NTHREADS]
                            [--gfaout]

UGMAGrefiner is a tool which 1) recruit unitigs to MAGs to improve
compeleteness of MAGs which assembled by metaSPAdes and binning by existing
binning tools, it is able to assign edges to multiple bins. 2) identify
whether there are multiple genomes mixed in one MAG and get the unitigs of the
different part of these genomes. UGMAGrefiner uses the connectivity and
coverage information from assembly graphs to improve completeness of MAGs from
exsiting binning tools on unitig level and to identify different part of
similar genomes mixed in one MAG.

optional arguments:
  -h, --help            show this help message and exit
  --edges EDGES         path to the assembly_graph.fastg file
  --graph GRAPH         path to the assembly_graph_after_simplification.gfa
                        file
  --paths PATHS         path to the contigs.paths file
  --binned BINNED       path to the .csv file with the initial binning output
                        from an existing tool
  --output OUTPUT       path to the output folder
  --prefix PREFIX       prefix for the output file
  --depth DEPTH         maximum depth for the breadth-first-search. [default:
                        5]
  --threshold THRESHOLD
                        threshold for determining inconsistent vertices.
                        [default: 1.5]
  --delimiter DELIMITER
                        delimiter for input/output results [default: ,
                        (comma)]
  --nthreads NTHREADS   number of threads to use. [default: 8]
  --gfaout              whether output gfa file for each MAG

# output
prefix_bin_group_Separated/New_cluster.fasta: the unitig sequences belong to "Separated" or "New"
prefix_bin*add.fasta: the newly added unitigs for each MAG
prefix_bin*add.gfa: the assembly graph in gfa format for each MAG with newly added unitigs
prefix_cluster.fasta: the sequences of all unitig clusters
prefix_extend_bin_cluster.csv: contain extended result for each unitig cluster 
prefix_newly_identified_cluster.csv: contain every unitigs' information of newly identified unitig clusters
prefix_UGMAGrefiner.log: log
prefix_UGMAGrefienr_output.csv: contain newly add unitigs for each MAG

