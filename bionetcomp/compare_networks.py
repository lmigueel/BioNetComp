import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pandas import DataFrame
import networkx as nx
import bionetcomp
import requests
import json

def enrich_string(genes,output_folder,name,fdr_cutoff,taxid):

	output = output_folder + '/'+name+'_enrichment.txt'
	report_file = open(output, "w")

	string_api_url = "https://string-db.org/api"	
	output_format = "json"
	method = "enrichment"


        ##
        ## Construct the request
        ##

	request_url = "/".join([string_api_url, output_format, method])

        ##
        ## Set parameters
        ##

	params = {

		"identifiers" : "%0d".join(genes), # your protein
		"species" : taxid, # species NCBI identifier 
		"caller_identity" : "www.awesome_app.org" # your app name

	}

        ##
        ## Call STRING
        ##

	response = requests.post(request_url, data=params)

        ##
        ## Read and parse the results
        ##

	data = json.loads(response.text)

	for row in data:

		term = row["term"]
		preferred_names = ",".join(row["preferredNames"])
		fdr = float(row["fdr"])
		description = row["description"]
		category = row["category"]

		if category == "Process" and fdr < fdr_cutoff:

                ## print significant GO Process annotations

			print("\t".join([term, preferred_names, str(fdr), description]),file=report_file)


def compare_networks(network1,network2,g1,g2,output_folder,fdr,taxid):

	######## 1. Create edges ################

	column_names = ["edge", "interaction_score"]
	network1_final = pd.DataFrame(columns = column_names)

	k=0
	for index, row in network1.iterrows():
	    #print(row['geneA'], row['geneB']
	    key_one = row['geneA'] + "--" + row['geneB']
	    key_two = row['geneB'] + "--" + row['geneA']	
	    network1_final.loc[k] = [key_one,row['interaction_score']]
	    k=k+1
	    network1_final.loc[k] = [key_two,row['interaction_score']]
	    k=k+1

	#network1_final	

	network2_final = pd.DataFrame(columns = column_names)
	
	k=0
	for index, row in network2.iterrows():
		key_one = row['geneA'] + "--" + row['geneB']
		key_two = row['geneB'] + "--" + row['geneA']
		network2_final.loc[k] = [key_one,row['interaction_score']]
		k=k+1
		network2_final.loc[k] = [key_two,row['interaction_score']]
		k=k+1
        #network2_final 

	############# 2. Nodes ########

	nodes1 = list(set(g1))
	nodes2 = list(set(g2))
	nodes_intersec =[]
	nodes_net2 = []
	for d in nodes2:
	    if d in nodes1:
	        nodes_intersec.append(d)
	    else:
	        nodes_net2.append(d)
        
	nodes_net1 = list(set(nodes1) - set(nodes_intersec))

	#print(len(nodes1),len(nodes2))
	#print(len(nodes_net1),len(nodes_net2),len(nodes_intersec))

	column_names = ["gene", "network1","network2","both","color"]
	node_report = pd.DataFrame(columns = column_names)

	all_genes = nodes_intersec+nodes_net1+nodes_net2
	#all_genes

	n1 = [1]*(len(nodes_intersec)+len(nodes_net1)) + [0]*(len(nodes_net2))
	n2 = [1]*len(nodes_intersec) + [0]*(len(nodes_net1)) + [1]*len(nodes_net2)
	b = [1]*(len(nodes_intersec)) + [0]*(len(nodes_net1)+len(nodes_net2))
	color_node = ['blue']*len(nodes_intersec) + ['red']*len(nodes_net1) + ['yellow']*len(nodes_net2)

	node_report['gene'] = all_genes
	node_report['network1'] = n1
	node_report['network2'] = n2
	node_report['both'] = b
	node_report['color'] = color_node

	output = output_folder + "/node_report.txt"
	node_report.to_csv(output,sep="\t",index=False)

	#enrichment - network1 exclusive nodes
	df_nodes1 = pd.DataFrame(nodes_net1) 
	output = output_folder + "/nodes_only_network1.txt"
	df_nodes1.to_csv(output,sep="\t",index=False,header=False)

	pass1=0
	pass2=0
	pass3=0

	if len(nodes_net1) > 0:	
		enrich_string(nodes_net1,output_folder,"enrichr_network1_exclusive_nodes",fdr,taxid)
		pass1=1
	else:
		print("Warning! No exclusive nodes in network1\n")

	#enrichment - network2 exclsuive nodes
	df_nodes2 = pd.DataFrame(nodes_net2)
	output = output_folder + "/nodes_only_network2.txt"
	df_nodes2.to_csv(output,sep="\t",index=False,header=False)

	if len(nodes_net2) > 0:
		enrich_string(nodes_net2,output_folder,"enrichr_network2_exclusive_nodes",fdr,taxid)
		pass2=1
	else:
		print("Warning! No exclusive nodes in network2\n")

	#enrichment - intersection nodes
	df_nodes12 = pd.DataFrame(nodes_intersec)
	output = output_folder + "/nodes_intersection_network12.txt"
	df_nodes12.to_csv(output,sep="\t",index=False,header=False)

	if len(nodes_intersec) > 0:
		enrich_string(nodes_intersec,output_folder,"enrichr_intersection_nodes",fdr,taxid)
		pass3=1
	else:
		print("Warning! No intersection nodes.\n")

	########### 3. Edges ###############

	
	edges1 = list(network1_final['edge'])
	edges2 = list(network2_final['edge'])

	edges_intersec =[]
	edges_net2 = []
	for e in edges2:
	    if e in edges1:
	        edges_intersec.append(e)
	    else:
	        edges_net2.append(e)

        
	edges_net1 = list(set(edges1) - set(edges_intersec))


	column_names = ["edge", "network1","network2","both","color"]
	edge_report = pd.DataFrame(columns = column_names)

	all_edges = edges_intersec+edges_net1+edges_net2

	e1 = [1]*(len(edges_intersec)+len(edges_net1)) + [0]*(len(edges_net2))
	e2 = [1]*len(edges_intersec) + [0]*(len(edges_net1)) + [1]*len(edges_net2)
	be = [1]*(len(edges_intersec)) + [0]*(len(edges_net1)+len(edges_net2))
	color_edge = ['red']*len(edges_intersec) + ['black']*len(edges_net1) + ['black']*len(edges_net2)

	edge_report['edge'] = all_edges
	edge_report['network1'] = e1
	edge_report['network2'] = e2
	edge_report['both'] = be
	edge_report['color'] = color_edge

	edge_report = edge_report.drop_duplicates()
	

	#generate networks from both networks
	#bionetcomp.edges_by_network(edge_report,output_folder,fdr,taxid)

	#### Network - edges1
	if pass1 ==1:
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
	else:
		print("Warning! No exclusive edges for network1.\n")


        #### Network - edges2

	if pass2 ==1:
		aux = edge_report.loc[(edge_report['network2']==1) & (edge_report['both']==0),]
		df_edges2 = pd.DataFrame(aux)
		aux1 = aux.loc[:,'edge'].str.split('--',n=1,expand=True)
		df_edges2['g1'] = aux1[0]
		df_edges2['g2'] = aux1[1]

		cols = df_edges2.columns.tolist()
		cols = ['g1','g2','edge', 'network1', 'network2', 'both', 'color']
		df_edges2 = df_edges2[cols]
		#df_edges2

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
	else:
		print("Warning! No exclusive edges for network2.\n")

	##### 4. Jaccard coefficient #########


        ## Nodes
	intersect_nodes = len(nodes_intersec)
	all_nodes = len(all_genes)
	jaccard_nodes = intersect_nodes/all_nodes

        ## edges
	intersect_edges = len(edges_intersec)/2
	all_edges2 = len(all_edges)/2
	jaccard_edges = intersect_edges/all_edges2

	labels = ["nodes","edges"]
	values = [jaccard_nodes,jaccard_edges]

	x = np.arange(len(labels))  # the label locations

	fig, ax = plt.subplots()
	rects1 = ax.bar(x, values, label='JC')

	def autolabel(rects):
		"""Attach a text label above each bar in *rects*, displaying its height."""
		for rect in rects:
			height = rect.get_height()
			ax.annotate('{}'.format(height),
			  xy=(rect.get_x() + rect.get_width() / 2, height),
			  xytext=(0, 3),  # 3 points vertical offset
			  textcoords="offset points",
			  ha='center', va='bottom')

        # Add some text for labels, title and custom x-axis tick labels, etc.
	ax.set_ylabel('Jaccard coefficient (JC)')
	ax.set_title('Jaccard coefficient plot - nodes and edges')	
	ax.set_xticks(x)
	ax.set_xticklabels(labels)
	ax.legend()

	autolabel(rects1)
	plt.ylim(0.0,max(values)+0.1)
	fig.tight_layout()

	output = output_folder + "/jaccard_coefficient_nodes_and_edges.png"
	plt.savefig(output)
	
	##### 5. final plot ########

	aux = edge_report
	df_final = pd.DataFrame(aux)
	aux2=aux.loc[:,'edge'].str.split('--',n=1,expand=True)
	df_final['g1'] = aux2[0]
	df_final['g2'] = aux2[1]

	cols = df_final.columns.tolist()
	cols = ['g1','g2','edge', 'network1', 'network2', 'both', 'color']
	df_final = df_final[cols]
	
	output = output_folder + "/edge_report.txt"
	df_final.to_csv(output,sep="\t",index=False)

	g= nx.from_pandas_edgelist(df_final,'g1','g2',edge_attr=True)

	#colors
	color_n = []
	for node in g:
	    c = node_report.loc[node_report['gene'] == node]['color'].values
	    color_n.append(c)
	    #print(n)
    
	color_e =[]
	for e in g.edges():
	    key = e[0]+"--"+e[1]
	    c = edge_report.loc[edge_report['edge'] == key]['color'].values
	    #print(c)
	    color_e.append(c)

	color_e = np.concatenate(color_e, axis=0)
	color_n = np.concatenate(color_n, axis=0)


	#FINAL PLOT
	pos = nx.spring_layout(g,k=0.30,iterations=20)
	betCent = nx.betweenness_centrality(g, normalized=True, endpoints=True)
	size =  [v * 10000 for v in betCent.values()]

	plt.figure(figsize=(20,20))

	nx.draw_networkx(g, pos=pos,with_labels=True,
        	         node_color=color_n,
                	 node_size=size,
	                 edge_color = np.squeeze(color_e))
	plt.axis('off')

	output = output_folder + "/final_network_plot.png"
	plt.savefig(output)	
