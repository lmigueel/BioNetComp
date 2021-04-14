import matplotlib
import matplotlib.pyplot as plt
import os
import pandas as pd
import numpy as np
from pandas import DataFrame
import networkx as nx
import collections

def comparative_metrics(network1,network2,output_folder):

	g1 = nx.from_pandas_edgelist(network1,'geneA','geneB',edge_attr=True)
	g2 = nx.from_pandas_edgelist(network2,'geneA','geneB',edge_attr=True)

	
	### 1. Nodes and edges ###
	nodes1 = g1.number_of_nodes()
	nodes2 = g2.number_of_nodes()
	edges1 = g1.number_of_edges()
	edges2 = g2.number_of_edges()

	labels = ["nodes","edges"]
	net1 = [nodes1,edges1]
	net2 = [nodes2,edges2]

	x = np.arange(len(labels))  # the label locations
	width = 0.35  # the width of the bars

	fig, ax = plt.subplots()
	rects1 = ax.bar(x - width/2, net1, width, label='Network1')
	rects2 = ax.bar(x + width/2, net2, width, label='Network2')

	# Add some text for labels, title and custom x-axis tick labels, etc.
	ax.set_ylabel('Count')
	ax.set_title('Basic metrics')
	ax.set_xticks(x)
	ax.set_xticklabels(labels)
	ax.legend()

	def autolabel(rects):
	    """Attach a text label above each bar in *rects*, displaying its height."""
	    for rect in rects:
	        height = rect.get_height()
        	ax.annotate('{}'.format(height),
                    xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom')

	autolabel(rects1)
	autolabel(rects2)

	fig.tight_layout()

	#plt.show()
	output = output_folder + "/nodes_edges.png"
	plt.savefig(output)

	##### 2. Closeness Centrality #######

	cc1 = nx.closeness_centrality(g1)
	cc2 = nx.closeness_centrality(g2)

	data = [list(cc1.values()), list(cc2.values())]
	fig1, ax1 = plt.subplots()
	ax1.set_title('Closeness Centrality boxplot')
	ax1.boxplot(data)
	plt.xticks([1, 2], ['network1', 'network2'])

	#plt.show()
	output = output_folder + "/closeness_centrality_boxplot.png"
	plt.savefig(output)

	cc1DF = pd.DataFrame(cc1.items(), columns=['Gene', 'Closeness Centrality'])
	cc1DF = cc1DF.sort_values('Closeness Centrality', ascending=False)
	output = output_folder + "/Closeness_Centrality_values_network1.txt"
	cc1DF.to_csv(output,sep="\t",index=False)

	cc2DF = pd.DataFrame(cc2.items(), columns=['Gene', 'Closeness Centrality'])
	cc2DF = cc2DF.sort_values('Closeness Centrality', ascending=False)
	output = output_folder + "/Closeness_Centrality_values_network2.txt"
	cc2DF.to_csv(output,sep="\t",index=False)	

	##### 3. Betweenness Centrality ######

	bc1 = nx.betweenness_centrality(g1)
	bc2 = nx.betweenness_centrality(g2)

	data = [list(bc1.values()), list(bc2.values())]
	fig2, ax2 = plt.subplots()
	ax2.set_title('Betweeness Centrality boxplot')
	ax2.boxplot(data)
	plt.xticks([1, 2], ['network1', 'network2'])

	#plt.show()
	output = output_folder + "/betweenness_centrality_boxplot.png"
	plt.savefig(output)

	bc1DF = pd.DataFrame(bc1.items(), columns=['Gene', 'Betweeness Centrality'])
	bc1DF = bc1DF.sort_values('Betweeness Centrality', ascending=False)
	output = output_folder + "/Betweeness_Centrality_values_network1.txt"
	bc1DF.to_csv(output,sep="\t",index=False)   

	bc2DF = pd.DataFrame(bc2.items(), columns=['Gene', 'Betweeness Centrality'])
	bc2DF = bc2DF.sort_values('Betweeness Centrality', ascending=False)
	output = output_folder + "/Betweeness_Centrality_values_network2.txt"
	bc2DF.to_csv(output,sep="\t",index=False)

	#### 4. Degree histogram ######

	degree_sequence = sorted([d for n, d in g1.degree()], reverse=True)  # degree sequence
	degreeCount = collections.Counter(degree_sequence)
	deg, cnt = zip(*degreeCount.items())

	plt.subplot(2,1,1)
	plt.bar(deg, cnt, width=0.80, color="b")
	plt.title("Degree Histogram - network1")
	plt.ylabel("Count")
	plt.xlabel("Degree")
	plt.xticks([d+0.0 for d in deg])

	plt.subplot(2,1,2)
	plt.title('Degree boxplot - network1')
	plt.boxplot(deg)
	plt.xticks([1], ['network1'])

	plt.tight_layout(pad=2.0)
	#plt.show()
	output = output_folder + "/network1_degree_hist.png"
	plt.savefig(output)

	plt.clf()
	degree_sequence = sorted([d for n, d in g2.degree()], reverse=True)  # degree sequence
	degreeCount = collections.Counter(degree_sequence)
	deg, cnt = zip(*degreeCount.items())

	plt.subplot(2,1,1)
	plt.bar(deg, cnt, width=0.80, color="b")
	plt.title("Degree Histogram - network2")
	plt.ylabel("Count")
	plt.xlabel("Degree")
	plt.xticks([d+0.0 for d in deg])

	plt.subplot(2,1,2)
	plt.title('Degree boxplot - network2')
	plt.boxplot(deg)
	plt.xticks([1], ['network2'])

	plt.tight_layout(pad=2.0)
        #plt.show()
	output = output_folder + "/network2_degree_hist.png"
	plt.savefig(output)

	######### 5. Metrics #############

	## reports txt
	output = output_folder + "/comparative_metrics_BioNetComp.txt"
	report_file = open(output, "w")

	print('### NETWORK 1  ###',file=report_file)
	print("Number of vertices:", g1.number_of_nodes(),file=report_file)
	print("Number of edges:", g1.number_of_edges(),file=report_file)
	print("Network is connected?",nx.is_connected(g1),file=report_file)

	if nx.is_connected(g1):
        	print('average path length:',nx.average_shortest_path_length(g1),file=report_file)
	        print('average diameter:',nx.diameter(g1),file=report_file)

	components = nx.connected_components(g1)
	largest_component = max(components, key=len)
	subgraph = g1.subgraph(largest_component)
	diameter = nx.diameter(subgraph)
	print("Network diameter of largest component:", diameter,file=report_file)

	print("Average clustering:",nx.average_clustering(g1),file=report_file)

	cc = nx.closeness_centrality(g1)
	df = pd.DataFrame.from_dict({
        	'node': list(cc.keys()),
	        'centrality': list(cc.values())
	    })
	selected = df.sort_values('centrality', ascending=False)['node'].iloc[0]
	value = df.sort_values('centrality', ascending=False)['centrality'].iloc[0]
	print("Node with max closeness centrality:[Node]-[Centrality]",selected,"-",value,file=report_file)

	bc = nx.betweenness_centrality(g1)
	df = pd.DataFrame.from_dict({
        	'node': list(bc.keys()),
	        'centrality': list(bc.values())
	    })
	selected = df.sort_values('centrality', ascending=False)['node'].iloc[0]
	value = df.sort_values('centrality', ascending=False)['centrality'].iloc[0]
	print("Node with max betweenness centrality:[Node]-[Centrality]",selected,"-",value,file=report_file)


	print("\n\n######### NETWORK2 ##########",file=report_file)
	print("Number of vertices:", g2.number_of_nodes(),file=report_file)
	print("Number of edges:", g2.number_of_edges(),file=report_file)
	print("Network is connected?",nx.is_connected(g2),file=report_file)

	if nx.is_connected(g2):
        	print('average path length:',nx.average_shortest_path_length(g2),file=report_file)
	        print('average diameter:',nx.diameter(g2),file=report_file)

	components = nx.connected_components(g2)
	largest_component = max(components, key=len)
	subgraph = g2.subgraph(largest_component)
	diameter = nx.diameter(subgraph)
	print("Network diameter of largest component:", diameter,file=report_file)

	print("Average clustering:",nx.average_clustering(g2),file=report_file)

	cc = nx.closeness_centrality(g2)
	df = pd.DataFrame.from_dict({
        	'node': list(cc.keys()),
	        'centrality': list(cc.values())
	    })
	selected = df.sort_values('centrality', ascending=False)['node'].iloc[0]
	value = df.sort_values('centrality', ascending=False)['centrality'].iloc[0]
	print("Node with max closeness centrality:[Node]-[Centrality]",selected,"-",value,file=report_file)

	bc = nx.betweenness_centrality(g2)
	df = pd.DataFrame.from_dict({
        	'node': list(bc.keys()),
	        'centrality': list(bc.values())
	    })
	selected = df.sort_values('centrality', ascending=False)['node'].iloc[0]
	value = df.sort_values('centrality', ascending=False)['centrality'].iloc[0]
	print("Node with max betweenness centrality:[Node]-[Centrality]",selected,"-",value,file=report_file)

	
	######  6. network plots ##########

	#network 1
	pos = nx.spring_layout(g1,k=0.30,iterations=20)
	betCent = nx.betweenness_centrality(g1, normalized=True, endpoints=True)
	color = [20000.0 * g1.degree(v) for v in g1]
	size =  [v * 10000 for v in betCent.values()]
	plt.figure(figsize=(20,20))
	nx.draw_networkx(g1, pos=pos, with_labels=True,
        	         node_color=color,
                	 node_size=size )
	plt.axis('off')
	output = output_folder + "/network1_plot.png"
	plt.savefig(output)

	#network 2
	pos = nx.spring_layout(g2,k=0.30,iterations=20)
	betCent = nx.betweenness_centrality(g2, normalized=True, endpoints=True)
	color = [20000.0 * g2.degree(v) for v in g2]
	size =  [v * 10000 for v in betCent.values()]
	plt.figure(figsize=(20,20))
	nx.draw_networkx(g2, pos=pos, with_labels=True,
                         node_color=color,
                         node_size=size )
	plt.axis('off')
	output = output_folder + "/network2_plot.png"
	plt.savefig(output)

	return g1,g2
