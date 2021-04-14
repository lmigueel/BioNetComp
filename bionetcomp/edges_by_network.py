import os
import pandas as pd
import numpy as np
from pandas import DataFrame
import networkx as nx
import requests 
import json
import matplotlib.pyplot as plt

def edges_by_network(edge_report,output_folder,fdr,taxid):

	#### Network - edges1
	aux = edge_report.loc[(edge_report['network1']==1) & (edge_report['both']==0),]
	df_edges1 = pd.DataFrame(aux)
	aux1 = aux.loc[:,'edge'].str.split('-',n=1,expand=True)
	df_edges1['g1'] = aux1[0]
	df_edges1['g2'] = aux1[1]

	cols = df_edges1.columns.tolist()
	cols = ['g1','g2','edge', 'network1', 'network2', 'both', 'color']
	df_edges1 = df_edges1[cols]
	#df_edges1

	ge1= nx.from_pandas_edgelist(df_edges1,'g1','g2',edge_attr=True)

	pos = nx.spring_layout(ge1,k=0.30,iterations=20)
	betCent = nx.betweenness_centrality(ge1, normalized=True, endpoints=True)
	size =  [v * 10000 for v in betCent.values()]

	plt.figure(figsize=(20,20))
	nx.draw_networkx(ge1, pos=pos,with_labels=True,
        	         node_color="blue",
                	 node_size=size,
	                 edge_color = "black")
	plt.axis('off')

	output = output_folder + "/only_network1_edges.png"
	plt.savefig(output)


	#### Network - edges2
	aux = edge_report.loc[(edge_report['network2']==1) & (edge_report['both']==0),]
	df_edges2 = pd.DataFrame(aux)
	aux1 = aux.loc[:,'edge'].str.split('--',n=1,expand=True)
	df_edges2['g1'] = aux1[0]
	df_edges2['g2'] = aux1[1]

	cols = df_edges2.columns.tolist()
	cols = ['g1','g2','edge', 'network1', 'network2', 'both', 'color']
	df_edges2 = df_edges2[cols]
        #df_edges1

	ge2= nx.from_pandas_edgelist(df_edges2,'g1','g2',edge_attr=True)

	pos = nx.spring_layout(ge2,k=0.30,iterations=20)
	betCent = nx.betweenness_centrality(ge2, normalized=True, endpoints=True)
	size =  [v * 10000 for v in betCent.values()]

	plt.figure(figsize=(20,20))
	nx.draw_networkx(ge2, pos=pos,with_labels=True,
                         node_color="yellow",
                         node_size=size,
                         edge_color = "black")
	plt.axis('off')

	output = output_folder + "/only_network2_edges.png"
	plt.savefig(output)
	
   	
