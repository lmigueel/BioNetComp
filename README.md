# BioNetComp

[![PyPI](https://img.shields.io/pypi/v/bionetcomp.svg?label=PyPI&color=green)](https://pypi.org/project/bionetcomp/)
[![Conda](https://img.shields.io/conda/vn/bioconda/bionetcomp.svg?label=Conda&color=green)](https://anaconda.org/bioconda/bionetcomp)
[![DOI](https://img.shields.io/badge/DOI-10.1101%2F2021.04.14.439897-red)](https://doi.org/10.1101/2021.04.14.439897)

BioNetComp: a Python package for biological network development and comparison 

        ██████╗ ██╗ ██████╗ ███╗   ██╗███████╗████████╗ ██████╗ ██████╗ ███╗   ███╗██████╗ 
        ██╔══██╗██║██╔═══██╗████╗  ██║██╔════╝╚══██╔══╝██╔════╝██╔═══██╗████╗ ████║██╔══██╗
        ██████╔╝██║██║   ██║██╔██╗ ██║█████╗     ██║   ██║     ██║   ██║██╔████╔██║██████╔╝
        ██╔══██╗██║██║   ██║██║╚██╗██║██╔══╝     ██║   ██║     ██║   ██║██║╚██╔╝██║██╔═══╝ 
        ██████╔╝██║╚██████╔╝██║ ╚████║███████╗   ██║   ╚██████╗╚██████╔╝██║ ╚═╝ ██║██║     
        ╚═════╝ ╚═╝ ╚═════╝ ╚═╝  ╚═══╝╚══════╝   ╚═╝    ╚═════╝ ╚═════╝ ╚═╝     ╚═╝╚═╝     
                                                                                   
        version 1.1




- [Installation](#installing)
- [Overview](#overview)
- [Workflow](#workflow)
- [Command-line interface](#command-line-interface)
- [Python library usage](#python-library-usage)
- [Examples](#examples)

# Installing

To build and install from source, run

```shell
python setup.py install
```
You can also install from pip with

```shell
pip install bionetcomp
``` 

# Overview

BioNetComp compares two biologial networks from a list of genes/proteins. It provides both an easy-to-use object-oriented Python API and a command-line interface (CLI) for network comparison and post-analysis. 

# Workflow

BioNetComp contains a flowchart designed to provide a structured comparative approach between two biological networks through the STRING database, as well as metrics and comparative reports, in addition to network visualizations. 

BioNetComp features include:

    1. A text file containing the list of nodes and total edges, differentiated by color and presence and absence on the network;
    2. A text file containing exclusive nodes of each network and those in common;
    3. Exclusive networks and a final network plot, containing comparative information;
    4. Network plot generated only by exclusive edges of each biological network;
    5. Comparative graphics of the number of nodes and edges;
    6. Exclusive comparison charts, such as the betweenness and closeness centrality boxplots;
    7. Degree histogram chart and its boxplot for each network;
    8. Enrichment of the entire network, but also exclusive and common nodes.
    9. Betweenness and closeness centrality gene ranking;
    10. Jaccard coefficient between networks applied to nodes and edges for dissimilarity observations.


# Command-line interface

BioNetComp can be executed from the command line using the bionetcomp command. It takes two list files, output folder name and organism taxid (based on STRING) as input and outputs a several analysis of the biological network generated. You also may change the STRING interaction score threshold and enrichment cutoff. 

```
usage: bionetcomp [-h] --in1 IN1 --in2 IN2 --taxid TAXID --output_folder
                  OUTPUT_FOLDER [--fdr FDR] [--threshold THRESHOLD]

optional arguments:
  -h, --help            show this help message and exit
  --in1 IN1             File containing a first gene list.
  --in2 IN2             File containing a second gene list.
  --taxid TAXID         STRING taxonomy ID. Ex: 9606
  --output_folder OUTPUT_FOLDER  Output folder
  --fdr FDR             FDR cutoff for enrichment analysis
  --threshold THRESHOLD Threshold for STRING interaction score

example: python3 bionetcomp.py --in1 list1.txt --in2 list2.txt --output output_folder --taxid 9606
```
In the output_folder, BioNetComp generates the following outputs:

Remember: 'network1' and 'network2' labels are associated to --in1 and --in2 list of genes, respectively. If there are no edges for network plots, or any nodes for enrichment, a warning will be shown for the user. 

    1. node_report.txt  all nodes of both networks divided by network
    2. edge_report.txt  all edges of both networks divided by network
    3. nodes_only_network1.txt  list of nodes present only in network1
    4. nodes_only_network2.txt  list of nodes present only in network2
    5. nodes_intersection_network12.txt  list of nodes present in both networks
    6. enrichr_network1_exclusive_nodes_enrichment.txt  nodes_only_network1.txt enrichment
    7. enrichr_network2_exclusive_nodes_enrichment.txt  nodes_only_network2.txt enrichment
    8. enrichr_intersection_nodes_enrichment.txt   intersection nodes enrichment
    9. enrichment_network1.txt  enrichment for all nodes in network1
    10. enrichment_network2.txt enrichment for all nodes in network2
    11. Betweeness_Centrality_values_network1.txt betweeness centrality value for all nodes of network1                          
    12. Betweeness_Centrality_values_network2.txt betweeness centrality value for all nodes of network2
    13. Closeness_Centrality_values_network1.txt  closeness centrality value for all nodes of network1 
    14. Closeness_Centrality_values_network2.txt  closeness centrality value for all nodes of network2
    15. comparative_metrics_BioNetComp.txt  basic network metrics for both networks  
    16. betweenness_centrality_boxplot.png  betweeness centrality boxplot
    17. closeness_centrality_boxplot.png  closeness centrality boxplot
    18. network1_plot.png   network1 plot
    19. network2_plot.png  network2 plot
    20. network1_degree_hist.png  degree histogram of network1
    21. network2_degree_hist.png  degree histogram of network2
    22. only_network1_edges.png   network developed from edges only present in network1
    23. only_network2_edges.png   network developed from edges only present in network2
    24. nodes_edges.png   barplot of comparative numbers of nodes and edges
    25. jaccard_coefficient_nodes_and_edges.png  jaccard coefficient plot for nodes and edges
    26. final_network_plot.png   final network plot colored by nodes and edges presence
        
Colors of the plots:

(i)   In the final network, nodes in red, blue, and yellow represent the intersection nodes, exclusive nodes of network1, and exclusive nodes of network2, respectively. 

(ii)  Edges in red represent the intersection ones. 

(iii) In basics network plots, the degree of each node delimits the size of the node in the final network plot, allowing visualization of the hub nodes. The viridis pallete was used. 
                             
# Python library usage

BioNetComp generates a folder in currently folder with output_folder name. 

To use as a Python library

```python

import bionetcomp
import os

# BioNetComp arguments
input1 = '/opt/data/input1.txt'
input2 = '/opt/data/input2.txt'
output_folder = 'bionetcomp_output'
taxid = 4932
threshold = 0.05
fdr = 0.05

#create the output folder
os.system("mkdir %s"%output_folder)

#develoment
network1 = bionetcomp.network_development(input1,taxid,threshold)
network2 = bionetcomp.network_development(input2,taxid,threshold)

#comparative metrics
g1, g2 = bionetcomp.comparative_metrics(network1,network2,output_folder)
        
#enrichment
bionetcomp.network_enrichment(input1,taxid,fdr,output_folder,"network1")
bionetcomp.network_enrichment(input2,taxid,fdr,output_folder,"network2")

#network comparison    
bionetcomp.compare_networks(network1,network2,g1,g2,output_folder,fdr,taxid)

```

# Examples

A main folder contains two lists with 83 and 75 genes of *S. cerevisiae* (in1.txt and in2.txt files)

```shell
bionetcomp --in1 in1.txt --in2 in2.txt --taxid 4932 --output_folder bionetcomp_teste
```



