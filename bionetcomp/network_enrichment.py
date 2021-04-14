import requests ## python -m pip install requests 
import json
import pandas as pd

def network_enrichment(input_file,taxid,fdr_cutoff,output_folder,name):

	genes_list = pd.read_csv(input_file, sep='\t',header=None)
	genes = genes_list.iloc[:,0]
	my_genes = genes.tolist()

	column_names = ["terms", "names","fdr","description"]
	enrichdf = pd.DataFrame(columns = column_names)
	
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

	    "identifiers" : "%0d".join(my_genes), # your protein
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
	k=0
	for row in data:

		term = row["term"]
		preferred_names = ",".join(row["preferredNames"])
		fdr = float(row["fdr"])
		description = row["description"]
		category = row["category"]

		if category == "Process" and fdr < fdr_cutoff:

			## print significant GO Process annotations
			enrichdf.loc[k]=[term, preferred_names, str(fdr), description]
			k=k+1

	output = output_folder + "/enrichment_" + name +".txt"
	enrichdf.to_csv(output,sep="\t",index=False,header=False)		
