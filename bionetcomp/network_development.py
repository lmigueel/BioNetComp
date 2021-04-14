import os
import pandas as pd
import numpy as np
from pandas import DataFrame
import networkx as nx
import requests

def network_development(input_file,taxid,threshold):

	genes_list = pd.read_csv(input_file, sep='\t',header=None)
	genes = genes_list.iloc[:,0]
	my_genes = genes.tolist()

	column_names = ["geneA", "geneB","interaction_score"]
	network = pd.DataFrame(columns = column_names)

	string_api_url = "https://string-db.org/api"
	output_format = "tsv-no-header"
	method = "network"

	##
	## Construct the request
	##

	request_url = "/".join([string_api_url, output_format, method])

	##
	## Set parameters
	##

	#print(my_genes)

	params = {

	    "identifiers" : "%0d".join(my_genes), # your protein
	    "species" : taxid, # species NCBI identifier 
	    "caller_identity" : "www.awesome_app.org", # your app name
	}

	#print(params)

	##
	## Call STRING
	##

	response = requests.post(request_url, data=params)

	##
	## Read and parse the results
	##

	k=0
	check_inter1={}
	for line in response.text.strip().split("\n"):
        
        	l = line.strip().split("\t")
	        p1, p2 = l[2], l[3]

	        ## filter the interaction according to experimental score
	        experimental_score = float(l[5])
	        key_one = p1 + "--" + p2
	        key_two = p2 + "--" + p1
        	if experimental_score >= threshold and (key_one not in check_inter1 and key_two not in check_inter1):
                
	                #print("\t".join([p1, p2, "experimentally confirmed (prob. %.3f)" % experimental_score]))
        	        network.loc[k] = [p1, p2, experimental_score]
                	k=k+1
	                check_inter1[key_one]=1
        	        check_inter1[key_two]=1


	return network
