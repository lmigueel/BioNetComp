# -*- coding: utf-8 -*-
#
#   This file is part of the tspex package, available at:
#   https://github.com/lucasmiguel/BioNetComp/
#
#   Tspex is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program. If not, see <https://www.gnu.org/licenses/>.
#
#   Contact: lucasmigueel@gmail.com

"""
Command-line interface for bionetcomp
"""
import os
import argparse
import sys
import bionetcomp


out = '''
	██████╗ ██╗ ██████╗ ███╗   ██╗███████╗████████╗ ██████╗ ██████╗ ███╗   ███╗██████╗ 
	██╔══██╗██║██╔═══██╗████╗  ██║██╔════╝╚══██╔══╝██╔════╝██╔═══██╗████╗ ████║██╔══██╗
	██████╔╝██║██║   ██║██╔██╗ ██║█████╗     ██║   ██║     ██║   ██║██╔████╔██║██████╔╝
	██╔══██╗██║██║   ██║██║╚██╗██║██╔══╝     ██║   ██║     ██║   ██║██║╚██╔╝██║██╔═══╝ 
	██████╔╝██║╚██████╔╝██║ ╚████║███████╗   ██║   ╚██████╗╚██████╔╝██║ ╚═╝ ██║██║     
	╚═════╝ ╚═╝ ╚═════╝ ╚═╝  ╚═══╝╚══════╝   ╚═╝    ╚═════╝ ╚═════╝ ╚═╝     ╚═╝╚═╝     
                                                                                   
	version 1.1
	'''

print(out)


def bionetcomp_cli(in1,in2, taxid, output_folder,fdr,threshold):

	""" Biological network development and analysis."""

	os.system("mkdir %s"%output_folder)

	
	print(". Step 1: Interactions downloading ...\n")
	
	network1 = bionetcomp.network_development(in1,taxid,threshold)
	network2 = bionetcomp.network_development(in2,taxid,threshold)
	
	print(".. Step 2: Generating comparative metrics ...\n")
	
	g1, g2 = bionetcomp.comparative_metrics(network1,network2,output_folder)
	

	print("... Step 3: Enrichment of networks is almost done! ...\n")
	
	bionetcomp.network_enrichment(in1,taxid,fdr,output_folder,"network1")
	bionetcomp.network_enrichment(in2,taxid,fdr,output_folder,"network2")

	print(".... Step 4: Comparing networks ...\n")
	
	bionetcomp.compare_networks(network1,network2,g1,g2,output_folder,fdr,taxid)
        
def main():


	usage_text = '''example:

	 python3 bionetcomp.py --in1 list1.txt --in2 list2.txt --output output_folder --taxid 9606

	'''

	my_parser = argparse.ArgumentParser(epilog=usage_text)
	my_parser.add_argument('--in1', type=str, help='File containing a first gene list.',required=True)
	my_parser.add_argument('--in2', type=str, help='File containing a second gene list.',required=True)
	my_parser.add_argument('--taxid', type=int,help='STRING taxonomy ID. Ex: 9606',required=True)
	my_parser.add_argument('--output_folder', type=str,help='Output folder',required=True)
	my_parser.add_argument('--fdr', type=float,default=0.05,help='FDR cutoff for enrichment analysis',required=False)
	my_parser.add_argument('--threshold', type=float,default=0.4,help='Threshold for STRING interaction score',required=False)

	args = my_parser.parse_args()    


	if len(sys.argv) < 4:
		parser.print_help()
		sys.exit(0)
	
	bionetcomp_cli(**vars(args))

